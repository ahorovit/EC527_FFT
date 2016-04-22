/*

	04/21/2016
	This version comes from https://equilibriumofnothing.wordpress.com/2013/10/14/algorithm-iterative-fft/


	g++ -O0 -std=c++0x -o FFT FFT_equil.cpp -lrt



*/



#include <climits>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <time.h> 
using namespace std;


int unityArray[2]; //stores unity step
//timing
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
void iterativeFFT(const vector<complex<double> >& primal, vector<complex<double> >& dual,const int P) 
{
  const int N = primal.size();
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
    const complex<double> unityRoot(cos(theta), sin(theta));

    // each higher level doubles the step size
    for (int offset = 0; offset < primal.size(); offset += unityArray[1]) {
      complex<double> omega = 1;

      // combine within a step segment (note only iterate over half step)
      for (int k = 0; k < unityArray[1]/2; k++) {
        const complex<double> u = dual[offset + k];

        const complex<double> t = omega * dual[offset + k + unityStep/2];
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

int main(int argc, char *argv[]) {
  struct timespec time1, time2, diffTime;
  struct timespec diff(struct timespec start, struct timespec end);  
  int clock_gettime(clockid_t clk_id, struct timespec *tp);

  // input number of coefficients
  //cout << "input number of coefficients" << endl;
  int N = 1048576*2*2*2*2;
  //cin >> N;

  // easy case - assume N is even power of 2
  const int P = lg(N);

  // check
  if (N != pown(P)) {
    cout << "error, " << N << " is not an even power of 2" << endl;
    exit(1);
  }

  // random coefficients for a polynomial (primal)
  vector<complex<double> > primal(N, 0);
  for (int i = 0; i < N; i++)
    primal[i] = i;

  // transformed dual
  vector<complex<double> > dual(N, 0);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
  iterativeFFT(primal, dual, P); // primal -> dual
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
  // print original primal coefficients
  //cout << "primal:" << endl;
  //for (int i = 0; i < N; i++)
   // cout << i << "\t" << primal[i] << endl;

  //cout << endl;

  // print dual coefficients
  //cout << "dul:" << endl;
 // for (int i = 0; i < N; i++)
  //  cout << i << "\t" << dual[i] << endl;

  // need another array for inverse
  vector<complex<double> > dualPrime(N, 0);

  // use -P as flag for inverse
  iterativeFFT(dual, dualPrime, -P); // dual -> primal

 // cout << endl;

  // print dual coefficients
  cout << "dualPrime:" << endl;
 // for (int i = 0; i < N; i++)
    cout << N << "\t" << dualPrime[N-1] << endl;

//Timing output
  diffTime = diff(time1, time2); 
  printf("%ld.%.9ld seconds\n", (long long)diffTime.tv_sec, diffTime.tv_nsec); 


  exit(0);
}
