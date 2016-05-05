/*
		EC527 Final Project --> FFT CUDA implementation
		Adin Horovitz, Monil Jhaveri, Evan Bowman

		Serial code is a slowed-down version of example found at:
			https://equilibriumofnothing.wordpress.com/2013/10/14/algorithm-iterative-fft/

		CUDA implementation uses example as a guide --> all STL functionality implemented using
		structs and direct complex value calculation :(

		to compile:

		nvcc -o FFT FFT_CUDA.cu

*/

#include <cstdio>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <climits>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <time.h> 
#include <stdint.h>
#include <thrust/complex.h>


#define N 1048576
#define MAX_THREAD 512
#define GIG 1000000000


typedef std::complex<double> cpx;
typedef std::vector<cpx> cVec;


/*		timing struct	*/
struct timespec diff(struct timespec start, struct timespec end)
{
  struct timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp;
}


/* 		-- Complex Struct for CUDA --		*/
typedef struct compDouble{
	double real;
	double imag;
} d_cpx;


/* Forward Declare Serial Functions */
uint32_t reverseBits(uint32_t i);
int lg(uint32_t i);
int pown(const int p);
int find_absP(const int P);
int is_inverse(bool inverse);

void iterativeFFT(const cVec & primal, cVec & dual,const int P);
void run_CUDA(d_cpx * d_in, d_cpx * d_out, int P, bool inverse);

void init_vec(cVec &data);
void init_dVec(d_cpx * data);
void zero_dVec(d_cpx * data);
void print_dVec(d_cpx * data);
void print_cVec(cVec &data);
void print_both(cVec& cpu, d_cpx * gpu);
int compare_out(cVec& cpu, d_cpx * gpu, double prec);


// Assertion to check for errors
#define CUDA_SAFE_CALL(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
	if (code != cudaSuccess) 
	{
		fprintf(stderr,"CUDA_SAFE_CALL: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}



/*		---		CUDA KERNELS		---		*/

__global__ void kernel_butterfly(d_cpx * input, d_cpx * output, int pwr2)
{
	// Identify indices to be swapped
	uint32_t thisIdx;
	uint32_t thatIdx;
//	int numBlocks = N/MAX_THREAD;

	if ( N <= MAX_THREAD)
		thisIdx = (uint32_t)threadIdx.x;
	else
		thisIdx = (MAX_THREAD * blockIdx.x) + threadIdx.x; 

	thatIdx = thisIdx;


	// SWAR to reverse bits of thisIdx
	register uint32_t mask = 0x55555555; // 0101...
	thatIdx = ((thatIdx & mask) << 1) | ((thatIdx >> 1) & mask);
	mask = 0x33333333; // 0011...
	thatIdx = ((thatIdx & mask) << 2) | ((thatIdx >> 2) & mask);
	mask = 0x0f0f0f0f; // 00001111...
	thatIdx = ((thatIdx & mask) << 4) | ((thatIdx >> 4) & mask);
	mask = 0x00ff00ff; // 0000000011111111...
	thatIdx = ((thatIdx & mask) << 8) | ((thatIdx >> 8) & mask);
	// 00000000000000001111111111111111 no need for mask
	thatIdx = (thatIdx << 16) | (thatIdx >> 16);

	thatIdx = thatIdx >> (32 - pwr2);

	//	printf("thisInd: %d		thatInd: %d\n", thisIdx, thatIdx);


	// Swap values
	output[thisIdx] = input[thatIdx];
	output[thatIdx] = input[thisIdx];

}


__global__ void kernel_FFTstage(d_cpx * input, int uStep, d_cpx uRoot, bool inverse)
{

	int numK = uStep / 2;

	// Allocate space in shared memory
//	__shared__ d_cpx[2 * MAX_THREAD];


	// Block position determines offset and k --> see IterativeFFT for parallel
	int offset = ((blockIdx.x * MAX_THREAD + threadIdx.x) / numK) * uStep;  
	int k = (blockIdx.x * MAX_THREAD + threadIdx.x) % numK;		

//	printf("[%d,%d] offset: %d, k: %d\n",threadIdx.x, threadIdx.y, offset, k);

	// omega = uRoot ^ k
	d_cpx omega, temp;
	omega.real = 1; omega.imag = 0;
	for(int i = 0; i < k; i++)
	{
		temp.real = (omega.real*uRoot.real) - (omega.imag*uRoot.imag);
		temp.imag = (omega.real*uRoot.imag) + (omega.imag*uRoot.real);
	
		omega = temp;
	}

//	printf("[%d,%d] uRoot: (%.3f,%.3f) omega: (%.3f,%.3f)\n", offset, k, uRoot.real, uRoot.imag, omega.real, omega.imag);


	// FFT stage change to input:

	int thisInd = offset + k;
	int thatInd = offset + k + uStep/2;

	temp = input[thatInd];
	d_cpx u = input[thisInd];

//	printf("thr[%d,%d] off=%d k=%d this:%d=(%.3f,%.3f) that:%d=(%.3f,%.3f)\n",blockIdx.x ,threadIdx.x, offset, k, 
//			thisInd, u.real, u.imag, thatInd, temp.real, temp.imag);



	// t = omega * temp
	d_cpx t;
	t.real = (omega.real*temp.real) - (omega.imag*temp.imag);
	t.imag = (omega.real*temp.imag) + (omega.imag*temp.real);

//	printf("[%d,%d] u: (%.3f,%.3f) t: (%.3f,%.3f)\n",offset, 
//			k, u.real, u.imag, t.real, t.imag);


	// Reuse temp for adding u + t
	temp.real = u.real + t.real;
	temp.imag = u.imag + t.imag;
	input[thisInd] = temp;

	// Reuse temp for subtracting u - t
	temp.real = u.real - t.real;
	temp.imag = u.imag - t.imag;
	input[thatInd] = temp;


	// If inverse FFT AND this is the last stage, divide through by N
	if (inverse && uStep == N/2)
	{
		input[thisInd].real /= N;
		input[thisInd].imag /= N;

		input[thatInd].real /= N;
		input[thatInd].imag /= N;
	}
}

/*		Globals		*/
int unityArray[2]; //stores unity step



/******************************************************************
**********                 BEGIN MAIN                   ***********
******************************************************************/
int main(int argc, char *argv[])
{

	printf("**************************************************\n");
	printf("*******     FFT CUDA Implementation       ********\n");
	printf("*******            N = %d                ********\n", N);
	printf("**************************************************\n");

	// Timing Structs
	struct timespec time1, time2, diffTime;
	struct timespec diff(struct timespec start, struct timespec end);  
	int clock_gettime(clockid_t clk_id, struct timespec *tp);
	cudaEvent_t start, stop;	
	float elapsed_gpu;

	// easy case - assume N is even power of 2
	const int P = lg(N);

	// check to be sure N is pwr of 2
	if (N != pown(P)) {
	std::cout << "error, " << N << " is not an even power of 2" << std::endl;
		exit(1);
	}


	/*	---		Begin Serial (slow) execution	--- */


	// Serial Vectors can use STL
	cVec ser_in(N, 0);
	cVec ser_out(N, 0);
	init_vec(ser_in); // initialize vector for serial implementation


	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
	iterativeFFT(ser_in, ser_out, P); 
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);


/*
	// need another array for inverse
	cVec dualPrime(N, 0);

	// use -P as flag for inverse
	iterativeFFT(ser_out, ser_in, -P); 

*/
	//Timing output
	diffTime = diff(time1, time2); 




	/* ---		 Begin CUDA FFT			 	--- */


	// kernel complex numbers/comp arrays must use structs
	CUDA_SAFE_CALL(cudaSetDevice(1));
	size_t allocSize = N * sizeof(d_cpx);
	d_cpx * h_in = (d_cpx *)malloc(allocSize);
	d_cpx * h_out = (d_cpx *)malloc(allocSize);

	// initialize device "vector"
	init_dVec(h_in);	
	zero_dVec(h_out);

	// Allocate space on Device
	d_cpx * d_in, *d_out;	// pointer to Device Memory
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_in, allocSize));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_out, allocSize));


	// Transfer input array to device
	CUDA_SAFE_CALL(cudaMemcpy(d_in, h_in, allocSize, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_out, h_out, allocSize, cudaMemcpyHostToDevice));


	// Start CUDA timer
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	run_CUDA(d_in, d_out, P, false);

	// Stop and destroy the timer
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsed_gpu, start, stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

/*
	// Perform Inverse FFT:  YOU MUST COMMENT MEMCOPY FROM D_OUT IF REVERSING!!!
  run_CUDA(d_out, d_in, P, true);
	CUDA_SAFE_CALL(cudaMemcpy(h_out, d_in, allocSize, cudaMemcpyDeviceToHost));
*/


	// Transfer the results back to the host
	CUDA_SAFE_CALL(cudaMemcpy(h_out, d_out, allocSize, cudaMemcpyDeviceToHost));

//	print_dVec(h_out);


	/************** Print Outputs *************/
//	printf("\tCPU\t\t\tGPU\n");

	/**** Free Memory ****/
	CUDA_SAFE_CALL(cudaFree(d_in));
	CUDA_SAFE_CALL(cudaFree(d_out));
	free(h_in);





	/*****  Print timing  ******/
	//      printf("%ld", (long int)((double)(CPG)*(double)
	//		 (GIG * time_stamp[j][i].tv_sec + time_stamp[j][i].tv_nsec)));

	double time_CPU = (double)((GIG * diffTime.tv_sec + diffTime.tv_nsec) / 10e6);

//	printf("\n\n\n  CPU time: %ld.%.9ld (msec)\n", (long long)(diffTime.tv_sec * 1000), diffTime.tv_nsec * 1000); 
//	printf("\n\n\n  CPU time: %.6f (msec)\n", time_CPU); 
	printf("\n\n\nCPU time(msec)\tGPU time(msec)\n"); 
	printf("%.5f\t%.5f\n", time_CPU, elapsed_gpu);
	printf("  #Errors: %d\n\n\n\n",   compare_out(ser_out, h_out, 0.0001));

	return 0;
}



/****************************************************************
						END MAIN
****************************************************************/


void run_CUDA(d_cpx * d_in, d_cpx * d_out, int P, bool inverse)
{

/*
#define N 16
#define MAX_THREAD 8
#define NUM_BLOCKS N/MAX_BLOCK
#define THREADS_PER_BLOCK N/NUM_BLOCKS
*/

	int num_blocks = N / MAX_THREAD;
	if (num_blocks < 1)
		num_blocks = 1;
				
	dim3 reordGrid(num_blocks);
	dim3 reordBlock(N / num_blocks);


  printf("N: %d, #blocks: %d, TPB: %d\n", N, num_blocks, N/num_blocks);

	// Call Reorder Kernel
	kernel_butterfly<<< reordGrid, reordBlock>>>(d_in, d_out, P);



	/* -- FFT_Stage only requires N/2 total threads -- */
	/* -- Call 1/2 as many blocks as butterfly kernel -- */

	int FFT_blocks;
	int threads_per_block;

	if (num_blocks <= 1)
	{
		FFT_blocks = 1;
		threads_per_block = N / 2; 
		if (threads_per_block > MAX_THREAD)
		{
			printf("threads_per_block cannot exceed %d\n", MAX_THREAD);
			exit(0);
		}
	}
	else
	{
		FFT_blocks = num_blocks / 2;
		threads_per_block = MAX_THREAD;		
	}

	dim3 FFT_grid(FFT_blocks);
	dim3 dimBlock(threads_per_block);

	printf("#blocks: %d, TPB: %d\n", FFT_blocks, threads_per_block);

	int uStep;
	double theta;
	d_cpx uRoot;


	// Execute FFT P (lg(N)) times
	for (int i = 1; i <= P; i++)
	{

		// Update stage parameters
		uStep = 0x1 << i;
		theta = (inverse ? -2 : 2) * M_PI / uStep;
		uRoot.real = cos(theta);
		uRoot.imag = sin(theta);


//		printf("\nbegin Stage %d\n", i);

//		printf("theta: %.3f\n", theta);

		// Call FFT stage kernel
		kernel_FFTstage<<< FFT_grid, dimBlock >>>(d_out, uStep , uRoot, inverse);


		// Stage must complete before next begins
		cudaThreadSynchronize();

//		printf("end FFT Stage %d\n", i);
	}
}



void print_cVec(cVec & data)
{
	for(int i = 0; i < N; i++)
		printf("[%d]: (%.3f, %.3f)\n", i, std::real(data[i]), std::imag(data[i]));
}


void print_dVec(d_cpx * data)
{
	for(int i = 0; i < N; i++)
		printf("[%d]: (%.3f, %.3f)\n", i, data[i].real, data[i].imag);
}

void print_both(cVec& cpu, d_cpx * gpu)
{

  for(int i = 0; i < N; i++)
  {
	  printf("[%d]: (%.3f, %.3f)\t\t(%.3f, %.3f)\n", i, real(cpu[i]), imag(cpu[i]),
				gpu[i].real, gpu[i].imag);
	}
}


int compare_out(cVec& cpu, d_cpx * gpu, double prec)
{
	int errors = 0;

  for(int i = 0; i < N; i++)
  {
		
		if( abs(real(cpu[i]) - gpu[i].real) > prec  || abs(imag(cpu[i]) - gpu[i].imag) > prec)
			errors++;
	}
	
	return errors;
}


// Initialize Vector for serial implementation
void init_vec(cVec & data)
{
	for (int i = 0; i < N; i++)
		data[i] = i;
}


// initialize Vector for CUDA implementation
void init_dVec(d_cpx * data)
{
	for (int i = 0; i < N; i++)
	{
		data[i].real = i;
		data[i].imag = 0;
	}
}

// Zero CUDA array
void zero_dVec(d_cpx * data)
{
	for (int i = 0; i < N; i++)
	{
		data[i].real = 0;
		data[i].imag = 0;
	}
}


/*	Reorder Input values in Butterfly order */
uint32_t reverseBits(uint32_t i) {
  register uint32_t mask = 0x55555555; // 0101...
  i = ((i & mask) << 1) | ((i >> 1) & mask);
  mask = 0x33333333; // 0011...
  i = ((i & mask) << 2) | ((i >> 2) & mask);
  mask = 0x0f0f0f0f; // 00001111...
  i = ((i & mask) << 4) | ((i >> 4) & mask);
  mask = 0x00ff00ff; // 0000000011111111...
  i = ((i & mask) << 8) | ((i >> 8) & mask);
  // 00000000000000001111111111111111 no need for mask
  i = (i << 16) | (i >> 16);
  return i;
}

int lg(uint32_t i) {
  int count = -1;
  while (i) {
    i = i >> 1;
    count++;
  }
  return count;
}

// Russian peasant algorithm
// Checks if input is pwr of two
int pown(const int p) {
  uint32_t w = p;
  w |= w >> 1;
  w |= w >> 2;
  w |= w >> 4;
  w |= w >> 8;
  w |= w >> 16;
  uint32_t mask = w & ~(w >> 1);

  int a = 1;
  while (mask) {
    a = a * a;
    if (mask & p)
      a *= 2;
    mask >>= 1;
  }

  return a;
}

//determines if the P is positive or negative
int find_absP(const int P){ 
	int temp = P; 
	if (temp < 0){
		return -temp; 
	} 	
	else return temp; 

}

 
int is_inverse(bool inverse){
	if (inverse){
		return -1;
	}
	else return 1;	
} 

// FFT takes complex input vector, container for complex output, and lg2(N) (will be negative if inverse)
void iterativeFFT(const cVec & primal, cVec & dual,const int P) 
{
//  const int N = primal.size();
  const bool inverse = P < 0;
  const int absP = find_absP(P);	// Whether positive or negative, absP is always positive

  // bottom level of iteration tree --> puts elements in butterfly order
  for (int i = 0; i < N; i++)
    dual[i] = primal[reverseBits(i) >> (32 - absP)];

  // there are absP levels above the bottom
  for (int p = 1; p <= find_absP(P); p++) {

    // complex root of unity
    const int unityStep = 0x1 << p;	// --> starts at two, doubles each iteration
    unityArray[1] = unityStep; 
    
 
    const double theta = (is_inverse(inverse)) * 2 * M_PI / unityArray[1]; // INVERSE
    const cpx unityRoot(cos(theta), sin(theta));

    // each higher level doubles the step size
    for (int offset = 0; offset < primal.size(); offset += unityArray[1]) {
      cpx omega = 1;

      // combine within a step segment (note only iterate over half step)
      for (int k = 0; k < unityArray[1]/2; k++) {
        cpx u = dual[offset + k];

        const cpx t = omega * dual[offset + k + unityStep/2];
        omega *= unityRoot;

        dual[offset + k] = u + t;
        dual[offset + k + unityArray[1]/2] = u - t;
      }
    }
  }

  if (inverse) // INVERSE
    for (int j = 0; j < primal.size(); j++)
      dual[j] /= N;
}

