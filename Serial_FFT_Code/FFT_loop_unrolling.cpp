/*

	04/21/2016
	This version comes from https://equilibriumofnothing.wordpress.com/2013/10/14/algorithm-iterative-fft/


	g++ -O0 -std=c++0x -o FFT FFT_basic_optimizations.cpp -lrt



*/



#include <climits>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <time.h> 
using namespace std;


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


inline uint32_t reverseBits(uint32_t i) {
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

void twiddle_factors(vector<complex<double> > & forwardTwiddle, vector<complex<double> > & backTwiddle, int N, const int P){

  const bool inverse = P < 0;
  const int absP = find_absP(P);	// Whether positive or negative, absP is always positive
  double twoPi = 2 * M_PI;  

  //double deltaFactor[absP]; // stores unityStep
  for (int i = 1; i <= absP; i++){
	
	const int unityStep = 0x1 << i;	// --> starts at two, doubles each iteration
	double theta =  twoPi/unityStep; 
	forwardTwiddle.push_back(complex<double>(cos(theta), sin(theta))); 
	backTwiddle.push_back(complex<double>(cos(-theta), sin(-theta)));  
   }
  int tempCounter = 0;  

  //const complex<double> unityRoot(cos(theta), sin(theta));
  //complex<double> twiddle1;    //for a different danileson switch cos and sin
} 	

inline void butterfly_loopunroll(const vector<complex<double> >& primal, vector<complex<double> >& dual, int N, int absP){ //unrolls the loop in the bit reversal based on the size of the input array 
	if (N > 32){
  		for ( int i = 0; i < N; i = i + 31){ 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+2] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+3] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+4] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+5] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+6] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+7] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+8] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+9] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+10] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+11] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+12] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+13] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+14] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+15] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+16] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+17] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+18] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+19] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+20] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+21] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+22] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+23] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+24] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+25] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+26] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+27] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+28] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+29] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+30] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+31] = primal[reverseBits(i+1) >> (32 - absP)];
  		} 
	} 
	else if (N > 16){
			
			int i = 0; 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+2] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+3] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+4] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+5] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+6] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+7] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+8] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+9] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+10] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+11] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+12] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+13] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+14] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+15] = primal[reverseBits(i+1) >> (32 - absP)];
		

	} 	
	else if (N > 8){ 

			int i = 0; 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+2] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+3] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+4] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+5] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+6] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+7] = primal[reverseBits(i+1) >> (32 - absP)];

	} 
	else if (N > 4){ 
			int i = 0; 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+2] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+3] = primal[reverseBits(i+1) >> (32 - absP)];
	} 	
	else{ 
			int i = 0; 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
	} 


} 	


inline void danielson_loopunroll(vector<complex<double> >& dual,  complex<double> unityRoot,  int unityStep, int offset){
      
	complex<double> omega = 1;

	if( unityStep < 32 ){

      		for (int k = 0; k < unityStep/2; k = k + 32) {
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep/2] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep/2 +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep/2 + 1] = u - t;	
			
        		 u = dual[offset + k + 2];
        		 t = omega * dual[offset + k + 2  + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 2] = u + t;
        		dual[offset + k + unityStep/2 + 2] = u - t;	
			
        		u = dual[offset + k + 3];
        		t = omega * dual[offset + k + 3  + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 3] = u + t;
        		dual[offset + k + unityStep/2 + 3] = u - t;	
			
        		u = dual[offset + k + 4];
        		t = omega * dual[offset + k + 4 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 4] = u + t;
        		dual[offset + k + unityStep/2 + 4] = u - t;	
			
        		u = dual[offset + k +  5];
        		t = omega * dual[offset + k + unityStep/2 + 5];
        		omega *= unityRoot;
        		dual[offset + k + 5] = u + t;
        		dual[offset + k + unityStep/2 + 5] = u - t;	
			
        		 u = dual[offset + k + 6];
        		 t = omega * dual[offset + k + 6  + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 6] = u + t;
        		dual[offset + k + unityStep/2 + 6] = u - t;	
			
        		u = dual[offset + k + 7];
        		t = omega * dual[offset + k + 7 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 7]  = u + t;
        		dual[offset + k + unityStep/2 + 7] = u - t;	
			
        		u = dual[offset + k + 8];
        		t = omega * dual[offset + k + unityStep/2 + 8];
        		omega *= unityRoot;
        		dual[offset + k + 8] = u + t;
        		dual[offset + k + unityStep/2 + 8] = u - t;	
			
        		u = dual[offset + k + 9];
        		t = omega * dual[offset + k + unityStep/2 + 9];
        		omega *= unityRoot;
        		dual[offset + k + 9] = u + t;
        		dual[offset + k + unityStep/2 + 9] = u - t;	
			
        		u = dual[offset + k + 10];
        		t = omega * dual[offset + k + unityStep/2 + 10];
        		omega *= unityRoot;
        		dual[offset + k + 10] = u + t;
        		dual[offset + k + unityStep/2 + 10] = u - t;	
			
        		u = dual[offset + k + 11];
        		t = omega * dual[offset + k +11 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 11] = u + t;
        		dual[offset + k + 11 + unityStep/2] = u - t;	
			
        		u = dual[offset + k + 12];
        		t = omega * dual[offset + k + 12 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k +12] = u + t;
        		dual[offset + k + unityStep/2 + 12] = u - t;	
			
        		u = dual[offset + k +13];
        		t = omega * dual[offset + k +13 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 13] = u + t;
        		dual[offset + k + unityStep/2 + 13] = u - t;	
			
        		u = dual[offset + k +14];
        		t = omega * dual[offset + k + unityStep/2 + 14];
        		omega *= unityRoot;
        		dual[offset + k + 14] = u + t;
        		dual[offset + k + 14 + unityStep/2] = u - t;	
			
        		u = dual[offset + k + 15];
        		t = omega * dual[offset + k + 15 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 15] = u + t;
        		dual[offset + k + unityStep/2 + 15] = u - t;	
			
        		u = dual[offset + k + 16];
        		t = omega * dual[offset + k + unityStep/2 + 16];
        		omega *= unityRoot;
        		dual[offset + k + 16] = u + t;
        		dual[offset + k + unityStep/2 + 16] = u - t;	
			
        		u = dual[offset + k + 17];
        		t = omega * dual[offset + k + unityStep/2 + 17];
        		omega *= unityRoot;
        		dual[offset + k + 17] = u + t;
        		dual[offset + k + unityStep/2 + 17] = u - t;	
			
        		u = dual[offset + k + 18 ];
        		t = omega * dual[offset + k +18 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 18] = u + t;
        		dual[offset + k + unityStep/2 + 18] = u - t;	
			
        		u = dual[offset + k + 19];
        		t = omega * dual[offset + k +19 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 19] = u + t;
        		dual[offset + k + unityStep/2 + 19] = u - t;	
			
        		u = dual[offset + k + 20];
        		t = omega * dual[offset + k + unityStep/2 + 20];
        		omega *= unityRoot;
        		dual[offset + k + 20] = u + t;
        		dual[offset + k + 20 + unityStep/2] = u - t;	
			
        		u = dual[offset + k + 21];
        		t = omega * dual[offset + k +21 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 21] = u + t;
        		dual[offset + k + unityStep/2 + 21] = u - t;	
			
        		u = dual[offset + k + 22];
        		t = omega * dual[offset + k + 22 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 22] = u + t;
        		dual[offset + k + unityStep/2 + 22] = u - t;	
			
        		u = dual[offset + k + 23];
        		t = omega * dual[offset + k + 23 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 23] = u + t;
        		dual[offset + k + unityStep/2 + 23] = u - t;	
			
        		u = dual[offset + k + 24];
        		t = omega * dual[offset + k + 24 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 24] = u + t;
        		dual[offset + k + unityStep/2 + 24 ] = u - t;	
			
        		u = dual[offset + k + 25];
        		t = omega * dual[offset + k + unityStep/2 + 25];
        		omega *= unityRoot;
        		dual[offset + k + 25 ] = u + t;
        		dual[offset + k + unityStep/2 + 25] = u - t;	
			
        		u = dual[offset + k + 26];
        		t = omega * dual[offset + k + 26 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 26] = u + t;
        		dual[offset + k + 26 + unityStep/2] = u - t;	
			
        		u = dual[offset + k + 27];
        		t = omega * dual[offset + k + 27 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 27 ] = u + t;
        		dual[offset + k + unityStep/2 + 27] = u - t;	
			
        		u = dual[offset + k + 28];
        		t = omega * dual[offset + k + unityStep/2 + 28];
        		omega *= unityRoot;
        		dual[offset + k + 28] = u + t;
        		dual[offset + k + unityStep/2 + 28] = u - t;	
			
        		u = dual[offset + k + 29];
        		t = omega * dual[offset + k + unityStep/2 + 29];
        		omega *= unityRoot;
        		dual[offset + k + 29] = u + t;
        		dual[offset + k + unityStep/2 + 29] = u - t;	
			
        		u = dual[offset + k + 30];
        		t = omega * dual[offset + k + unityStep/2 + 30];
        		omega *= unityRoot;
        		dual[offset + k + 30] = u + t;
        		dual[offset + k + unityStep/2 + 30] = u - t;	
			
        		u = dual[offset + k + 31];
        		t = omega * dual[offset + k + 31 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 31] = u + t;
        		dual[offset + k + unityStep/2 + 31] = u - t;
		} 
	} 
	else if( unityStep > 16){

			int k = 0; 
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep/2] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep/2 +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep/2 + 1] = u - t;	
			
        		 u = dual[offset + k + 2];
        		 t = omega * dual[offset + k + 2  + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 2] = u + t;
        		dual[offset + k + unityStep/2 + 2] = u - t;	
			
        		u = dual[offset + k + 3];
        		t = omega * dual[offset + k + 3  + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 3] = u + t;
        		dual[offset + k + unityStep/2 + 3] = u - t;	
			
        		u = dual[offset + k + 4];
        		t = omega * dual[offset + k + 4 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 4] = u + t;
        		dual[offset + k + unityStep/2 + 4] = u - t;	
			
        		u = dual[offset + k +  5];
        		t = omega * dual[offset + k + unityStep/2 + 5];
        		omega *= unityRoot;
        		dual[offset + k + 5] = u + t;
        		dual[offset + k + unityStep/2 + 5] = u - t;	
			
        		 u = dual[offset + k + 6];
        		 t = omega * dual[offset + k + 6  + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 6] = u + t;
        		dual[offset + k + unityStep/2 + 6] = u - t;	
			
        		u = dual[offset + k + 7];
        		t = omega * dual[offset + k + 7 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 7]  = u + t;
        		dual[offset + k + unityStep/2 + 7] = u - t;	
			
        		u = dual[offset + k + 8];
        		t = omega * dual[offset + k + unityStep/2 + 8];
        		omega *= unityRoot;
        		dual[offset + k + 8] = u + t;
        		dual[offset + k + unityStep/2 + 8] = u - t;	
			
        		u = dual[offset + k + 9];
        		t = omega * dual[offset + k + unityStep/2 + 9];
        		omega *= unityRoot;
        		dual[offset + k + 9] = u + t;
        		dual[offset + k + unityStep/2 + 9] = u - t;	
			
        		u = dual[offset + k + 10];
        		t = omega * dual[offset + k + unityStep/2 + 10];
        		omega *= unityRoot;
        		dual[offset + k + 10] = u + t;
        		dual[offset + k + unityStep/2 + 10] = u - t;	
			
        		u = dual[offset + k + 11];
        		t = omega * dual[offset + k +11 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 11] = u + t;
        		dual[offset + k + 11 + unityStep/2] = u - t;	
			
        		u = dual[offset + k + 12];
        		t = omega * dual[offset + k + 12 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k +12] = u + t;
        		dual[offset + k + unityStep/2 + 12] = u - t;	
			
        		u = dual[offset + k +13];
        		t = omega * dual[offset + k +13 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 13] = u + t;
        		dual[offset + k + unityStep/2 + 13] = u - t;	
			
        		u = dual[offset + k +14];
        		t = omega * dual[offset + k + unityStep/2 + 14];
        		omega *= unityRoot;
        		dual[offset + k + 14] = u + t;
        		dual[offset + k + 14 + unityStep/2] = u - t;	
			
        		u = dual[offset + k + 15];
        		t = omega * dual[offset + k + 15 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 15] = u + t;
        		dual[offset + k + unityStep/2 + 15] = u - t;	

	} 
	else if (unityStep > 7){ 

			int k = 0; 
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep/2] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep/2 +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep/2 + 1] = u - t;	
			
        		 u = dual[offset + k + 2];
        		 t = omega * dual[offset + k + 2  + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 2] = u + t;
        		dual[offset + k + unityStep/2 + 2] = u - t;	
			
        		u = dual[offset + k + 3];
        		t = omega * dual[offset + k + 3  + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 3] = u + t;
        		dual[offset + k + unityStep/2 + 3] = u - t;	
			
        		u = dual[offset + k + 4];
        		t = omega * dual[offset + k + 4 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 4] = u + t;
        		dual[offset + k + unityStep/2 + 4] = u - t;	
			
        		u = dual[offset + k +  5];
        		t = omega * dual[offset + k + unityStep/2 + 5];
        		omega *= unityRoot;
        		dual[offset + k + 5] = u + t;
        		dual[offset + k + unityStep/2 + 5] = u - t;	
			
        		 u = dual[offset + k + 6];
        		 t = omega * dual[offset + k + 6  + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 6] = u + t;
        		dual[offset + k + unityStep/2 + 6] = u - t;	
			
        		u = dual[offset + k + 7];
        		t = omega * dual[offset + k + 7 + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 7]  = u + t;
        		dual[offset + k + unityStep/2 + 7] = u - t;	
	} 
	else if ( unityStep > 3){ 
			int k = 0; 
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep/2] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep/2 +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep/2 + 1] = u - t;	
			
        		 u = dual[offset + k + 2];
        		 t = omega * dual[offset + k + 2  + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 2] = u + t;
        		dual[offset + k + unityStep/2 + 2] = u - t;	
			
        		u = dual[offset + k + 3];
        		t = omega * dual[offset + k + 3  + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k + 3] = u + t;
        		dual[offset + k + unityStep/2 + 3] = u - t;	


	} 
	else{ 

			int k = 0; 
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep/2];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep/2] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep/2 +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep/2 + 1] = u - t;	

	} 
} 



void loopunrollFFT(const vector<complex<double> >& primal, vector<complex<double> >& dual, const vector<complex<double> >& directionTwiddle  ,const int P) 
{
  const int N = primal.size();
  const bool inverse = P < 0;
  
  const int absP = find_absP(P);	// Whether positive or negative, absP is always positive
  butterfly_loopunroll(primal, dual, N, absP); 
  // bottom level of iteration tree --> puts elements in butterfly order
  /*for (int i = 0; i < N; i++){ 
    dual[i] = primal[reverseBits(i) >> (32 - absP)];
    dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
    dual[i+2] = primal[reverseBits(i+2) >> (32 - absP)];
    dual[i+3] = primal[reverseBits(i+3) >> (32 - absP)];
  } 
  */
  //start of the danielson 
  // there are absP levels above the bottom
  for (int p = 1; p <= absP; p++) {

    // complex root of unity
    const int unityStep = 0x1 << p;	// --> starts at two, doubles each iteration
    
 
    //const double theta = (inverse ? -1 : 1) * 2 * M_PI / unityStep; // INVERSE
    complex<double> unityRoot = directionTwiddle[p-1] ;

    // each higher level doubles the step size
    for (int offset = 0; offset < N; offset += unityStep) {

      // combine within a step segment (note only iterate over half step)
 	danielson_loopunroll(dual, unityRoot, unityStep/2, offset); 
/*      for (int k = 0; k < unityStep/2; k++) {
        const complex<double> u = dual[offset + k];

        const complex<double> t = omega * dual[offset + k + unityStep/2];
        omega *= unityRoot;

        dual[offset + k] = u + t;
        dual[offset + k + unityStep/2] = u - t;
      }
      */
    }
   
  }

  if (inverse) // INVERSE
    for (int j = 0; j < primal.size(); j++)
      dual[j] /= N;
}


void twiddleFFT(const vector<complex<double> >& primal, vector<complex<double> >& dual, const vector<complex<double> >& directionTwiddle  ,const int P) 
{
  const int N = primal.size();
  const bool inverse = P < 0;
  
  const int absP = find_absP(P);	// Whether positive or negative, absP is always positive

  // bottom level of iteration tree --> puts elements in butterfly order
  for (int i = 0; i < N; i++)
    dual[i] = primal[reverseBits(i) >> (32 - absP)];

  //start of the danielson 
  // there are absP levels above the bottom
  for (int p = 1; p <= absP; p++) {

    // complex root of unity
    const int unityStep = 0x1 << p;	// --> starts at two, doubles each iteration
    
 
    //const double theta = (inverse ? -1 : 1) * 2 * M_PI / unityStep; // INVERSE
    complex<double> unityRoot = directionTwiddle[p-1] ;

    // each higher level doubles the step size
    for (int offset = 0; offset < N; offset += unityStep) {
      complex<double> omega = 1;

      // combine within a step segment (note only iterate over half step)
      for (int k = 0; k < unityStep/2; k++) {
        const complex<double> u = dual[offset + k];

        const complex<double> t = omega * dual[offset + k + unityStep/2];
        omega *= unityRoot;

        dual[offset + k] = u + t;
        dual[offset + k + unityStep/2] = u - t;
      }
    }
  }

  if (inverse) // INVERSE
    for (int j = 0; j < primal.size(); j++)
      dual[j] /= N;
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
  for (int p = 1; p <= absP; p++) {

    // complex root of unity
    const int unityStep = 0x1 << p;	// --> starts at two, doubles each iteration
    
 
    const double theta = (inverse ? -1 : 1) * 2 * M_PI / unityStep; // INVERSE
    const complex<double> unityRoot(cos(theta), sin(theta));

    // each higher level doubles the step size
    for (int offset = 0; offset < N; offset += unityStep) {
      complex<double> omega = 1;

      // combine within a step segment (note only iterate over half step)
      for (int k = 0; k < unityStep/2; k++) {
        const complex<double> u = dual[offset + k];

        const complex<double> t = omega * dual[offset + k + unityStep/2];
        omega *= unityRoot;

        dual[offset + k] = u + t;
        dual[offset + k + unityStep/2] = u - t;
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
  
//  int N = 1048576*2*2*2*2; //use this for actual results 
 //int N = 16; 
  int N = 512; 
 
  //cin >> N;

  // easy case - assume N is even power of 2
  const int P = lg(N);
  vector<complex<double> > forwardTwiddle;
  vector<complex<double> > backTwiddle;
  twiddle_factors(forwardTwiddle, backTwiddle, N, P); 
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
  twiddle_factors(forwardTwiddle, backTwiddle, N, P); 
  twiddleFFT(primal, dual, forwardTwiddle, P); // primal -> dual
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
