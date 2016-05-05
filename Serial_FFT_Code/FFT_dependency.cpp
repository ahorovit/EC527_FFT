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
#define GIG 1000000000
#define CPG 1.596           // Cycles per GHz -- Adjust to your computer

using namespace std;

int unityArray[2];  

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


void baselineFFT(const vector<complex<double> >& primal, vector<complex<double> >& dual,const int P) 
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


inline void butterfly_loopunroll(const vector<complex<double> >& primal, vector<complex<double> >& dual, int N, int absP){ //unrolls the loop in the bit reversal based on the size of the input array 
	if (N > 31){
  		for ( int i = 0; i < N; i = i + 32){ 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+2] = primal[reverseBits(i+2) >> (32 - absP)];
    			dual[i+3] = primal[reverseBits(i+3) >> (32 - absP)];
    			dual[i+4] = primal[reverseBits(i+4) >> (32 - absP)];
    			dual[i+5] = primal[reverseBits(i+5) >> (32 - absP)];
    			dual[i+6] = primal[reverseBits(i+6) >> (32 - absP)];
    			dual[i+7] = primal[reverseBits(i+7) >> (32 - absP)];
    			dual[i+8] = primal[reverseBits(i+8) >> (32 - absP)];
    			dual[i+9] = primal[reverseBits(i+9) >> (32 - absP)];
    			dual[i+10] = primal[reverseBits(i+10) >> (32 - absP)];
    			dual[i+11] = primal[reverseBits(i+11) >> (32 - absP)];
    			dual[i+12] = primal[reverseBits(i+12) >> (32 - absP)];
    			dual[i+13] = primal[reverseBits(i+13) >> (32 - absP)];
    			dual[i+14] = primal[reverseBits(i+14) >> (32 - absP)];
    			dual[i+15] = primal[reverseBits(i+15) >> (32 - absP)];
    			dual[i+16] = primal[reverseBits(i+16) >> (32 - absP)];
    			dual[i+17] = primal[reverseBits(i+17) >> (32 - absP)];
    			dual[i+18] = primal[reverseBits(i+18) >> (32 - absP)];
    			dual[i+19] = primal[reverseBits(i+19) >> (32 - absP)];
    			dual[i+20] = primal[reverseBits(i+20) >> (32 - absP)];
    			dual[i+21] = primal[reverseBits(i+21) >> (32 - absP)];
    			dual[i+22] = primal[reverseBits(i+22) >> (32 - absP)];
    			dual[i+23] = primal[reverseBits(i+23) >> (32 - absP)];
    			dual[i+24] = primal[reverseBits(i+24) >> (32 - absP)];
    			dual[i+25] = primal[reverseBits(i+25) >> (32 - absP)];
    			dual[i+26] = primal[reverseBits(i+26) >> (32 - absP)];
    			dual[i+27] = primal[reverseBits(i+27) >> (32 - absP)];
    			dual[i+28] = primal[reverseBits(i+28) >> (32 - absP)];
    			dual[i+29] = primal[reverseBits(i+29) >> (32 - absP)];
    			dual[i+30] = primal[reverseBits(i+30) >> (32 - absP)];
    			dual[i+31] = primal[reverseBits(i+31) >> (32 - absP)];
  		} 
	} 
	else if (N > 16){
			
			int i = 0; 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+2] = primal[reverseBits(i+2) >> (32 - absP)];
    			dual[i+3] = primal[reverseBits(i+3) >> (32 - absP)];
    			dual[i+4] = primal[reverseBits(i+4) >> (32 - absP)];
    			dual[i+5] = primal[reverseBits(i+5) >> (32 - absP)];
    			dual[i+6] = primal[reverseBits(i+6) >> (32 - absP)];
    			dual[i+7] = primal[reverseBits(i+7) >> (32 - absP)];
    			dual[i+8] = primal[reverseBits(i+8) >> (32 - absP)];
    			dual[i+9] = primal[reverseBits(i+9) >> (32 - absP)];
    			dual[i+10] = primal[reverseBits(i+10) >> (32 - absP)];
    			dual[i+11] = primal[reverseBits(i+11) >> (32 - absP)];
    			dual[i+12] = primal[reverseBits(i+12) >> (32 - absP)];
    			dual[i+13] = primal[reverseBits(i+13) >> (32 - absP)];
    			dual[i+14] = primal[reverseBits(i+14) >> (32 - absP)];
    			dual[i+15] = primal[reverseBits(i+15) >> (32 - absP)];
		

	} 	
	else if (N > 8){ 

			int i = 0; 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+2] = primal[reverseBits(i+2) >> (32 - absP)];
    			dual[i+3] = primal[reverseBits(i+3) >> (32 - absP)];
    			dual[i+4] = primal[reverseBits(i+4) >> (32 - absP)];
    			dual[i+5] = primal[reverseBits(i+5) >> (32 - absP)];
    			dual[i+6] = primal[reverseBits(i+6) >> (32 - absP)];
    			dual[i+7] = primal[reverseBits(i+7) >> (32 - absP)];

	} 
	else if (N > 4){ 
			int i = 0; 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+2] = primal[reverseBits(i+2) >> (32 - absP)];
    			dual[i+3] = primal[reverseBits(i+3) >> (32 - absP)];
	} 	
	else{ 
			int i = 0; 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
	} 


} 	

inline void butterfly_loopunroll16(const vector<complex<double> >& primal, vector<complex<double> >& dual, int N, int absP){ //unrolls the loop in the bit reversal based on the size of the input array 
	if (N > 15){
			
  		for ( int i = 0; i < N; i = i + 16){ 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+2] = primal[reverseBits(i+2) >> (32 - absP)];
    			dual[i+3] = primal[reverseBits(i+3) >> (32 - absP)];
    			dual[i+4] = primal[reverseBits(i+4) >> (32 - absP)];
    			dual[i+5] = primal[reverseBits(i+5) >> (32 - absP)];
    			dual[i+6] = primal[reverseBits(i+6) >> (32 - absP)];
    			dual[i+7] = primal[reverseBits(i+7) >> (32 - absP)];
    			dual[i+8] = primal[reverseBits(i+8) >> (32 - absP)];
    			dual[i+9] = primal[reverseBits(i+9) >> (32 - absP)];
    			dual[i+10] = primal[reverseBits(i+10) >> (32 - absP)];
    			dual[i+11] = primal[reverseBits(i+11) >> (32 - absP)];
    			dual[i+12] = primal[reverseBits(i+12) >> (32 - absP)];
    			dual[i+13] = primal[reverseBits(i+13) >> (32 - absP)];
    			dual[i+14] = primal[reverseBits(i+14) >> (32 - absP)];
    			dual[i+15] = primal[reverseBits(i+15) >> (32 - absP)];
		}
		

	} 	
	else if (N > 7){ 

			int i = 0; 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+2] = primal[reverseBits(i+2) >> (32 - absP)];
    			dual[i+3] = primal[reverseBits(i+3) >> (32 - absP)];
    			dual[i+4] = primal[reverseBits(i+4) >> (32 - absP)];
    			dual[i+5] = primal[reverseBits(i+5) >> (32 - absP)];
    			dual[i+6] = primal[reverseBits(i+6) >> (32 - absP)];
    			dual[i+7] = primal[reverseBits(i+7) >> (32 - absP)];

	} 
	else if (N > 3){ 
			int i = 0; 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+2] = primal[reverseBits(i+2) >> (32 - absP)];
    			dual[i+3] = primal[reverseBits(i+3) >> (32 - absP)];
	} 	
	else{ 
			int i = 0; 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
	} 


} 	
inline void butterfly_loopunroll8(const vector<complex<double> >& primal, vector<complex<double> >& dual, int N, int absP){ //unrolls the loop in the bit reversal based on the size of the input array 
	if (N > 7){ 

  		for ( int i = 0; i < N; i = i + 8){ 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+2] = primal[reverseBits(i+2) >> (32 - absP)];
    			dual[i+3] = primal[reverseBits(i+3) >> (32 - absP)];
    			dual[i+4] = primal[reverseBits(i+4) >> (32 - absP)];
    			dual[i+5] = primal[reverseBits(i+5) >> (32 - absP)];
    			dual[i+6] = primal[reverseBits(i+6) >> (32 - absP)];
    			dual[i+7] = primal[reverseBits(i+7) >> (32 - absP)];

		} 
	} 
	else if (N > 3){ 
			int i = 0; 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+2] = primal[reverseBits(i+2) >> (32 - absP)];
    			dual[i+3] = primal[reverseBits(i+3) >> (32 - absP)];
	} 	
	else{ 
			int i = 0; 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
	} 


} 	
inline void butterfly_loopunroll4(const vector<complex<double> >& primal, vector<complex<double> >& dual, int N, int absP){ //unrolls the loop in the bit reversal based on the size of the input array 
	if (N > 3){ 
  		for ( int i = 0; i < N; i = i + 4){ 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
    			dual[i+2] = primal[reverseBits(i+2) >> (32 - absP)];
    			dual[i+3] = primal[reverseBits(i+3) >> (32 - absP)];
		}
	} 	
	else{ 
			int i = 0; 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
	} 


} 	
inline void butterfly_loopunroll2(const vector<complex<double> >& primal, vector<complex<double> >& dual, int N, int absP){ //unrolls the loop in the bit reversal based on the size of the input array 
  		for ( int i = 0; i < N; i = i + 2){ 
			dual[i] = primal[reverseBits(i) >> (32 - absP)];
    			dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
	} 


} 	

inline void danielson_loopunroll(vector<complex<double> >& dual,  complex<double> unityRoot,  int unityStep, int offset){
      
	complex<double> omega = 1;

	if( unityStep >  31 ){
      		for (int k = 0; k < unityStep; k = k + 32) {
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep + 1] = u - t;	
			
        		 u = dual[offset + k + 2];
        		 t = omega * dual[offset + k + 2  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 2] = u + t;
        		dual[offset + k + unityStep + 2] = u - t;	
			
        		u = dual[offset + k + 3];
        		t = omega * dual[offset + k + 3  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 3] = u + t;
        		dual[offset + k + unityStep + 3] = u - t;	
			
        		u = dual[offset + k + 4];
        		t = omega * dual[offset + k + 4 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 4] = u + t;
        		dual[offset + k + unityStep + 4] = u - t;	
			
        		u = dual[offset + k +  5];
        		t = omega * dual[offset + k + unityStep + 5];
        		omega *= unityRoot;
        		dual[offset + k + 5] = u + t;
        		dual[offset + k + unityStep + 5] = u - t;	
			
        		 u = dual[offset + k + 6];
        		 t = omega * dual[offset + k + 6  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 6] = u + t;
        		dual[offset + k + unityStep + 6] = u - t;	
			
        		u = dual[offset + k + 7];
        		t = omega * dual[offset + k + 7 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 7]  = u + t;
        		dual[offset + k + unityStep + 7] = u - t;	
			
        		u = dual[offset + k + 8];
        		t = omega * dual[offset + k + unityStep + 8];
        		omega *= unityRoot;
        		dual[offset + k + 8] = u + t;
        		dual[offset + k + unityStep + 8] = u - t;	
			
        		u = dual[offset + k + 9];
        		t = omega * dual[offset + k + unityStep + 9];
        		omega *= unityRoot;
        		dual[offset + k + 9] = u + t;
        		dual[offset + k + unityStep + 9] = u - t;	
			
        		u = dual[offset + k + 10];
        		t = omega * dual[offset + k + unityStep + 10];
        		omega *= unityRoot;
        		dual[offset + k + 10] = u + t;
        		dual[offset + k + unityStep + 10] = u - t;	
			
        		u = dual[offset + k + 11];
        		t = omega * dual[offset + k +11 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 11] = u + t;
        		dual[offset + k + 11 + unityStep] = u - t;	
			
        		u = dual[offset + k + 12];
        		t = omega * dual[offset + k + 12 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k +12] = u + t;
        		dual[offset + k + unityStep + 12] = u - t;	
			
        		u = dual[offset + k +13];
        		t = omega * dual[offset + k +13 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 13] = u + t;
        		dual[offset + k + unityStep + 13] = u - t;	
			
        		u = dual[offset + k +14];
        		t = omega * dual[offset + k + unityStep + 14];
        		omega *= unityRoot;
        		dual[offset + k + 14] = u + t;
        		dual[offset + k + 14 + unityStep] = u - t;	
			
        		u = dual[offset + k + 15];
        		t = omega * dual[offset + k + 15 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 15] = u + t;
        		dual[offset + k + unityStep + 15] = u - t;	
			
        		u = dual[offset + k + 16];
        		t = omega * dual[offset + k + unityStep + 16];
        		omega *= unityRoot;
        		dual[offset + k + 16] = u + t;
        		dual[offset + k + unityStep + 16] = u - t;	
			
        		u = dual[offset + k + 17];
        		t = omega * dual[offset + k + unityStep + 17];
        		omega *= unityRoot;
        		dual[offset + k + 17] = u + t;
        		dual[offset + k + unityStep + 17] = u - t;	
			
        		u = dual[offset + k + 18 ];
        		t = omega * dual[offset + k +18 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 18] = u + t;
        		dual[offset + k + unityStep + 18] = u - t;	
			
        		u = dual[offset + k + 19];
        		t = omega * dual[offset + k +19 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 19] = u + t;
        		dual[offset + k + unityStep + 19] = u - t;	
			
        		u = dual[offset + k + 20];
        		t = omega * dual[offset + k + unityStep + 20];
        		omega *= unityRoot;
        		dual[offset + k + 20] = u + t;
        		dual[offset + k + 20 + unityStep] = u - t;	
			
        		u = dual[offset + k + 21];
        		t = omega * dual[offset + k +21 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 21] = u + t;
        		dual[offset + k + unityStep + 21] = u - t;	
			
        		u = dual[offset + k + 22];
        		t = omega * dual[offset + k + 22 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 22] = u + t;
        		dual[offset + k + unityStep + 22] = u - t;	
			
        		u = dual[offset + k + 23];
        		t = omega * dual[offset + k + 23 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 23] = u + t;
        		dual[offset + k + unityStep + 23] = u - t;	
			
        		u = dual[offset + k + 24];
        		t = omega * dual[offset + k + 24 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 24] = u + t;
        		dual[offset + k + unityStep + 24 ] = u - t;	
			
        		u = dual[offset + k + 25];
        		t = omega * dual[offset + k + unityStep + 25];
        		omega *= unityRoot;
        		dual[offset + k + 25 ] = u + t;
        		dual[offset + k + unityStep + 25] = u - t;	
			
        		u = dual[offset + k + 26];
        		t = omega * dual[offset + k + 26 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 26] = u + t;
        		dual[offset + k + 26 + unityStep] = u - t;	
			
        		u = dual[offset + k + 27];
        		t = omega * dual[offset + k + 27 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 27 ] = u + t;
        		dual[offset + k + unityStep + 27] = u - t;	
			
        		u = dual[offset + k + 28];
        		t = omega * dual[offset + k + unityStep + 28];
        		omega *= unityRoot;
        		dual[offset + k + 28] = u + t;
        		dual[offset + k + unityStep + 28] = u - t;	
			
        		u = dual[offset + k + 29];
        		t = omega * dual[offset + k + unityStep + 29];
        		omega *= unityRoot;
        		dual[offset + k + 29] = u + t;
        		dual[offset + k + unityStep + 29] = u - t;	
			
        		u = dual[offset + k + 30];
        		t = omega * dual[offset + k + unityStep + 30];
        		omega *= unityRoot;
        		dual[offset + k + 30] = u + t;
        		dual[offset + k + unityStep + 30] = u - t;	
			
        		u = dual[offset + k + 31];
        		t = omega * dual[offset + k + 31 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 31] = u + t;
        		dual[offset + k + unityStep + 31] = u - t;
		} 
	} 

	else if( unityStep > 15){
 //     		for (int k = 0; k < unityStep; k = k + 16) {
			int k = 0; 
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep + 1] = u - t;	
			
        		 u = dual[offset + k + 2];
        		 t = omega * dual[offset + k + 2  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 2] = u + t;
        		dual[offset + k + unityStep + 2] = u - t;	
			
        		u = dual[offset + k + 3];
        		t = omega * dual[offset + k + 3  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 3] = u + t;
        		dual[offset + k + unityStep + 3] = u - t;	
			
        		u = dual[offset + k + 4];
        		t = omega * dual[offset + k + 4 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 4] = u + t;
        		dual[offset + k + unityStep + 4] = u - t;	
			
        		u = dual[offset + k +  5];
        		t = omega * dual[offset + k + unityStep + 5];
        		omega *= unityRoot;
        		dual[offset + k + 5] = u + t;
        		dual[offset + k + unityStep + 5] = u - t;	
			
        		 u = dual[offset + k + 6];
        		 t = omega * dual[offset + k + 6  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 6] = u + t;
        		dual[offset + k + unityStep + 6] = u - t;	
			
        		u = dual[offset + k + 7];
        		t = omega * dual[offset + k + 7 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 7]  = u + t;
        		dual[offset + k + unityStep + 7] = u - t;	
			
        		u = dual[offset + k + 8];
        		t = omega * dual[offset + k + unityStep + 8];
        		omega *= unityRoot;
        		dual[offset + k + 8] = u + t;
        		dual[offset + k + unityStep + 8] = u - t;	
			
        		u = dual[offset + k + 9];
        		t = omega * dual[offset + k + unityStep + 9];
        		omega *= unityRoot;
        		dual[offset + k + 9] = u + t;
        		dual[offset + k + unityStep + 9] = u - t;	
			
        		u = dual[offset + k + 10];
        		t = omega * dual[offset + k + unityStep + 10];
        		omega *= unityRoot;
        		dual[offset + k + 10] = u + t;
        		dual[offset + k + unityStep + 10] = u - t;	
			
        		u = dual[offset + k + 11];
        		t = omega * dual[offset + k +11 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 11] = u + t;
        		dual[offset + k + 11 + unityStep] = u - t;	
			
        		u = dual[offset + k + 12];
        		t = omega * dual[offset + k + 12 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k +12] = u + t;
        		dual[offset + k + unityStep + 12] = u - t;	
			
        		u = dual[offset + k +13];
        		t = omega * dual[offset + k +13 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 13] = u + t;
        		dual[offset + k + unityStep + 13] = u - t;	
			
        		u = dual[offset + k +14];
        		t = omega * dual[offset + k + unityStep + 14];
        		omega *= unityRoot;
        		dual[offset + k + 14] = u + t;
        		dual[offset + k + 14 + unityStep] = u - t;	
			
        		u = dual[offset + k + 15];
        		t = omega * dual[offset + k + 15 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 15] = u + t;
        		dual[offset + k + unityStep + 15] = u - t;	
		

	} 
	else if (unityStep > 7){ 
		
			int k = 0; 
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep + 1] = u - t;	
			
        		 u = dual[offset + k + 2];
        		 t = omega * dual[offset + k + 2  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 2] = u + t;
        		dual[offset + k + unityStep + 2] = u - t;	
			
        		u = dual[offset + k + 3];
        		t = omega * dual[offset + k + 3  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 3] = u + t;
        		dual[offset + k + unityStep + 3] = u - t;	
			
        		u = dual[offset + k + 4];
        		t = omega * dual[offset + k + 4 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 4] = u + t;
        		dual[offset + k + unityStep + 4] = u - t;	
			
        		u = dual[offset + k +  5];
        		t = omega * dual[offset + k + unityStep + 5];
        		omega *= unityRoot;
        		dual[offset + k + 5] = u + t;
        		dual[offset + k + unityStep + 5] = u - t;	
			
        		 u = dual[offset + k + 6];
        		 t = omega * dual[offset + k + 6  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 6] = u + t;
        		dual[offset + k + unityStep + 6] = u - t;	
			
        		u = dual[offset + k + 7];
        		t = omega * dual[offset + k + 7 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 7]  = u + t;
        		dual[offset + k + unityStep + 7] = u - t;	
		} 
	 
	else if ( unityStep > 3){ 
	
			int k = 0; 
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep + 1] = u - t;	
			
        		 u = dual[offset + k + 2];
        		 t = omega * dual[offset + k + 2  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 2] = u + t;
        		dual[offset + k + unityStep + 2] = u - t;	
			
        		u = dual[offset + k + 3];
        		t = omega * dual[offset + k + 3  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 3] = u + t;
        		dual[offset + k + unityStep + 3] = u - t;	
		 

	} 
	else if (unityStep > 1){ 
			int k = 0; 
			

        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep + 1] = u - t;	

	} 
	else{ 
			int k = 0; 
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
	} 

} 


inline void danielson_loopunroll16(vector<complex<double> >& dual,  complex<double> unityRoot,  int unityStep, int offset){
      
	complex<double> omega = 1;

	 if( unityStep > 15){
     		for (int k = 0; k < unityStep; k = k + 16) {
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep + 1] = u - t;	
			
        		 u = dual[offset + k + 2];
        		 t = omega * dual[offset + k + 2  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 2] = u + t;
        		dual[offset + k + unityStep + 2] = u - t;	
			
        		u = dual[offset + k + 3];
        		t = omega * dual[offset + k + 3  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 3] = u + t;
        		dual[offset + k + unityStep + 3] = u - t;	
			
        		u = dual[offset + k + 4];
        		t = omega * dual[offset + k + 4 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 4] = u + t;
        		dual[offset + k + unityStep + 4] = u - t;	
			
        		u = dual[offset + k +  5];
        		t = omega * dual[offset + k + unityStep + 5];
        		omega *= unityRoot;
        		dual[offset + k + 5] = u + t;
        		dual[offset + k + unityStep + 5] = u - t;	
			
        		 u = dual[offset + k + 6];
        		 t = omega * dual[offset + k + 6  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 6] = u + t;
        		dual[offset + k + unityStep + 6] = u - t;	
			
        		u = dual[offset + k + 7];
        		t = omega * dual[offset + k + 7 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 7]  = u + t;
        		dual[offset + k + unityStep + 7] = u - t;	
			
        		u = dual[offset + k + 8];
        		t = omega * dual[offset + k + unityStep + 8];
        		omega *= unityRoot;
        		dual[offset + k + 8] = u + t;
        		dual[offset + k + unityStep + 8] = u - t;	
			
        		u = dual[offset + k + 9];
        		t = omega * dual[offset + k + unityStep + 9];
        		omega *= unityRoot;
        		dual[offset + k + 9] = u + t;
        		dual[offset + k + unityStep + 9] = u - t;	
			
        		u = dual[offset + k + 10];
        		t = omega * dual[offset + k + unityStep + 10];
        		omega *= unityRoot;
        		dual[offset + k + 10] = u + t;
        		dual[offset + k + unityStep + 10] = u - t;	
			
        		u = dual[offset + k + 11];
        		t = omega * dual[offset + k +11 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 11] = u + t;
        		dual[offset + k + 11 + unityStep] = u - t;	
			
        		u = dual[offset + k + 12];
        		t = omega * dual[offset + k + 12 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k +12] = u + t;
        		dual[offset + k + unityStep + 12] = u - t;	
			
        		u = dual[offset + k +13];
        		t = omega * dual[offset + k +13 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 13] = u + t;
        		dual[offset + k + unityStep + 13] = u - t;	
			
        		u = dual[offset + k +14];
        		t = omega * dual[offset + k + unityStep + 14];
        		omega *= unityRoot;
        		dual[offset + k + 14] = u + t;
        		dual[offset + k + 14 + unityStep] = u - t;	
			
        		u = dual[offset + k + 15];
        		t = omega * dual[offset + k + 15 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 15] = u + t;
        		dual[offset + k + unityStep + 15] = u - t;	
		} 
		

	} 
	else if (unityStep > 7){ 
		
			int k = 0; 
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep + 1] = u - t;	
			
        		 u = dual[offset + k + 2];
        		 t = omega * dual[offset + k + 2  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 2] = u + t;
        		dual[offset + k + unityStep + 2] = u - t;	
			
        		u = dual[offset + k + 3];
        		t = omega * dual[offset + k + 3  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 3] = u + t;
        		dual[offset + k + unityStep + 3] = u - t;	
			
        		u = dual[offset + k + 4];
        		t = omega * dual[offset + k + 4 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 4] = u + t;
        		dual[offset + k + unityStep + 4] = u - t;	
			
        		u = dual[offset + k +  5];
        		t = omega * dual[offset + k + unityStep + 5];
        		omega *= unityRoot;
        		dual[offset + k + 5] = u + t;
        		dual[offset + k + unityStep + 5] = u - t;	
			
        		 u = dual[offset + k + 6];
        		 t = omega * dual[offset + k + 6  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 6] = u + t;
        		dual[offset + k + unityStep + 6] = u - t;	
			
        		u = dual[offset + k + 7];
        		t = omega * dual[offset + k + 7 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 7]  = u + t;
        		dual[offset + k + unityStep + 7] = u - t;	
		} 
	 
	else if ( unityStep > 3){ 
	
			int k = 0; 
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep + 1] = u - t;	
			
        		 u = dual[offset + k + 2];
        		 t = omega * dual[offset + k + 2  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 2] = u + t;
        		dual[offset + k + unityStep + 2] = u - t;	
			
        		u = dual[offset + k + 3];
        		t = omega * dual[offset + k + 3  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 3] = u + t;
        		dual[offset + k + unityStep + 3] = u - t;	
		 

	} 
	else if (unityStep > 1){ 
			int k = 0; 
			

        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep + 1] = u - t;	

	} 
	else{ 
			int k = 0; 
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
	} 

} 

inline void danielson_loopunroll8(vector<complex<double> >& dual,  complex<double> unityRoot,  int unityStep, int offset){
      
	complex<double> omega = 1;

	if (unityStep > 7){ 
		
     		for (int k = 0; k < unityStep; k = k + 8) {
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep + 1] = u - t;	
			
        		 u = dual[offset + k + 2];
        		 t = omega * dual[offset + k + 2  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 2] = u + t;
        		dual[offset + k + unityStep + 2] = u - t;	
			
        		u = dual[offset + k + 3];
        		t = omega * dual[offset + k + 3  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 3] = u + t;
        		dual[offset + k + unityStep + 3] = u - t;	
			
        		u = dual[offset + k + 4];
        		t = omega * dual[offset + k + 4 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 4] = u + t;
        		dual[offset + k + unityStep + 4] = u - t;	
			
        		u = dual[offset + k +  5];
        		t = omega * dual[offset + k + unityStep + 5];
        		omega *= unityRoot;
        		dual[offset + k + 5] = u + t;
        		dual[offset + k + unityStep + 5] = u - t;	
			
        		 u = dual[offset + k + 6];
        		 t = omega * dual[offset + k + 6  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 6] = u + t;
        		dual[offset + k + unityStep + 6] = u - t;	
			
        		u = dual[offset + k + 7];
        		t = omega * dual[offset + k + 7 + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 7]  = u + t;
        		dual[offset + k + unityStep + 7] = u - t;	
		} 
	} 
	 
	else if ( unityStep > 3){ 
	
			int k = 0; 
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep + 1] = u - t;	
			
        		 u = dual[offset + k + 2];
        		 t = omega * dual[offset + k + 2  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 2] = u + t;
        		dual[offset + k + unityStep + 2] = u - t;	
			
        		u = dual[offset + k + 3];
        		t = omega * dual[offset + k + 3  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 3] = u + t;
        		dual[offset + k + unityStep + 3] = u - t;	
		 

	} 
	else if (unityStep > 1){ 
			int k = 0; 
			

        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep + 1] = u - t;	

	} 
	else{ 
			int k = 0; 
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
	} 

} 
inline void danielson_loopunroll4(vector<complex<double> >& dual,  complex<double> unityRoot,  int unityStep, int offset){
      
	complex<double> omega = 1;
	
	if ( unityStep > 3){ 
	
     		for (int k = 0; k < unityStep; k = k + 4) {
			int k = 0; 
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep + 1] = u - t;	
			
        		 u = dual[offset + k + 2];
        		 t = omega * dual[offset + k + 2  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 2] = u + t;
        		dual[offset + k + unityStep + 2] = u - t;	
			
        		u = dual[offset + k + 3];
        		t = omega * dual[offset + k + 3  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 3] = u + t;
        		dual[offset + k + unityStep + 3] = u - t;	
		} 
		 

	} 
	else if (unityStep > 1){ 
			int k = 0; 
			

        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep + 1] = u - t;	

	} 
	else{ 
			int k = 0; 
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
	} 

} 
inline void danielson_loopunroll2(vector<complex<double> >& dual,  complex<double> unityRoot,  int unityStep, int offset){
      
	complex<double> omega = 1;
	
	if (unityStep > 1){ 
     		for (int k = 0; k < unityStep; k = k + 2) {
			

        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep + 1] = u - t;	

		} 
	} 
	else{ 
			int k = 0; 
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
	} 

} 

inline void danielson_loopunroll_dependency(vector<complex<double> >& dual,  complex<double> unityRoot,  int unityStep, int offset){
      
	complex<double> omega = 1;
	complex<double> omega1 = 1;
	complex<double> omega2 = 1;
	complex<double> omega3 = 1;
	complex<double> omega4 = 1;
	complex<double> omega5 = 1;
	complex<double> omega6 = 1;
	complex<double> omega7 = 1;
	complex<double> omega8 = 1;

	if (unityStep > 7){ 
		
     		for (int k = 0; k < unityStep; k = k + 8) {
			
			omega1 = unityRoot;  
			omega2 = pow(unityRoot, 2);
                        omega3 = pow(unityRoot, 3); 
			omega4 =  pow(unityRoot, 4);
                        omega5 = pow(unityRoot, 5);
                        omega6 = pow(unityRoot, 6);
                        omega7 = pow(unityRoot, 7);
                        omega8 = pow(unityRoot, 8);

        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega1 * dual[offset + k + unityStep +1];
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep + 1] = u - t;	
			
        		 u = dual[offset + k + 2];
        		 t = omega2 * dual[offset + k + 2  + unityStep];
        		dual[offset + k + 2] = u + t;
        		dual[offset + k + unityStep + 2] = u - t;	
			
        		u = dual[offset + k + 3];
        		t = omega3 * dual[offset + k + 3  + unityStep];
        		dual[offset + k + 3] = u + t;
        		dual[offset + k + unityStep + 3] = u - t;	
			
        		u = dual[offset + k + 4];
        		t = omega4 * dual[offset + k + 4 + unityStep];
        		dual[offset + k + 4] = u + t;
        		dual[offset + k + unityStep + 4] = u - t;	
			
        		u = dual[offset + k +  5];
        		t = omega5 * dual[offset + k + unityStep + 5];
        		omega *= unityRoot;
        		dual[offset + k + 5] = u + t;
        		dual[offset + k + unityStep + 5] = u - t;	
			
        		 u = dual[offset + k + 6];
        		 t = omega6 * dual[offset + k + 6  + unityStep];
        		dual[offset + k + 6] = u + t;
        		dual[offset + k + unityStep + 6] = u - t;	
			
        		u = dual[offset + k + 7];
        		t = omega7 * dual[offset + k + 7 + unityStep];
        		dual[offset + k + 7]  = u + t;
        		dual[offset + k + unityStep + 7] = u - t;	
		} 
	} 
	 
	else if ( unityStep > 3){ 
	
			int k = 0; 
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep + 1] = u - t;	
			
        		 u = dual[offset + k + 2];
        		 t = omega * dual[offset + k + 2  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 2] = u + t;
        		dual[offset + k + unityStep + 2] = u - t;	
			
        		u = dual[offset + k + 3];
        		t = omega * dual[offset + k + 3  + unityStep];
        		omega *= unityRoot;
        		dual[offset + k + 3] = u + t;
        		dual[offset + k + unityStep + 3] = u - t;	
		 

	} 
	else if (unityStep > 1){ 
			int k = 0; 
			

        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
			

        		u = dual[offset + k + 1];
        		t = omega * dual[offset + k + unityStep +1];
        		omega *= unityRoot;
        		dual[offset + k + 1] = u + t;
        		dual[offset + k + unityStep + 1] = u - t;	

	} 
	else{ 
			int k = 0; 
        		complex<double> u = dual[offset + k];
        		complex<double> t = omega * dual[offset + k + unityStep];
        		omega *= unityRoot;
        		dual[offset + k] = u + t;
        		dual[offset + k + unityStep] = u - t;	
	} 

} 
void loopunrollFFT_dep(const vector<complex<double> >& primal, vector<complex<double> >& dual, const vector<complex<double> >& directionTwiddle  ,const int P) 
{
  const int N = primal.size();
  const bool inverse = P < 0;
  
  const int absP = find_absP(P);	// Whether positive or negative, absP is always positive
  butterfly_loopunroll(primal, dual, N, absP); 
  // bottom level of iteration tree --> puts elements in butterfly order
/* for (int i = 0; i < N; i++){ 
    dual[i] = primal[reverseBits(i) >> (32 - absP)];
  //  dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
   // dual[i+2] = primal[reverseBits(i+2) >> (32 - absP)];
    //dual[i+3] = primal[reverseBits(i+3) >> (32 - absP)];
  } 
 */ 
  
  //start of the danielson 
  // there are absP levels above the bottom
  for (int p = 1; p <= absP; p++) {

    // complex root of unity
    const int unityStep = 0x1 << p;	// --> starts at two, doubles each iteration
    
   // complex<double> omega = 1;
 
    //const double theta = (inverse ? -1 : 1) * 2 * M_PI / unityStep; // INVERSE
    complex<double> unityRoot = directionTwiddle[p-1] ;

    // each higher level doubles the step size
    for (int offset = 0; offset < N; offset += unityStep) {

   //    combine within a step segment (note only iterate over half step)
   danielson_loopunroll_dependency(dual, directionTwiddle[p-1], unityStep/2, offset); 
  /*    for (int k = 0; k < unityStep/2; k++) {
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

void loopunrollFFT_32(const vector<complex<double> >& primal, vector<complex<double> >& dual, const vector<complex<double> >& directionTwiddle  ,const int P) 
{
  const int N = primal.size();
  const bool inverse = P < 0;
  
  const int absP = find_absP(P);	// Whether positive or negative, absP is always positive
  butterfly_loopunroll(primal, dual, N, absP); 
  // bottom level of iteration tree --> puts elements in butterfly order
/* for (int i = 0; i < N; i++){ 
    dual[i] = primal[reverseBits(i) >> (32 - absP)];
  //  dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
   // dual[i+2] = primal[reverseBits(i+2) >> (32 - absP)];
    //dual[i+3] = primal[reverseBits(i+3) >> (32 - absP)];
  } 
 */ 
  
  //start of the danielson 
  // there are absP levels above the bottom
  for (int p = 1; p <= absP; p++) {

    // complex root of unity
    const int unityStep = 0x1 << p;	// --> starts at two, doubles each iteration
    
   // complex<double> omega = 1;
 
    //const double theta = (inverse ? -1 : 1) * 2 * M_PI / unityStep; // INVERSE
    complex<double> unityRoot = directionTwiddle[p-1] ;

    // each higher level doubles the step size
    for (int offset = 0; offset < N; offset += unityStep) {

   //    combine within a step segment (note only iterate over half step)
   danielson_loopunroll(dual, directionTwiddle[p-1], unityStep/2, offset); 
  /*    for (int k = 0; k < unityStep/2; k++) {
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

void loopunrollFFT_16(const vector<complex<double> >& primal, vector<complex<double> >& dual, const vector<complex<double> >& directionTwiddle  ,const int P) 
{
  const int N = primal.size();
  const bool inverse = P < 0;
  
  const int absP = find_absP(P);	// Whether positive or negative, absP is always positive
  butterfly_loopunroll16(primal, dual, N, absP); 
  // bottom level of iteration tree --> puts elements in butterfly order
/* for (int i = 0; i < N; i++){ 
    dual[i] = primal[reverseBits(i) >> (32 - absP)];
  //  dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
   // dual[i+2] = primal[reverseBits(i+2) >> (32 - absP)];
    //dual[i+3] = primal[reverseBits(i+3) >> (32 - absP)];
  } 
 */ 
  
  //start of the danielson 
  // there are absP levels above the bottom
  for (int p = 1; p <= absP; p++) {

    // complex root of unity
    const int unityStep = 0x1 << p;	// --> starts at two, doubles each iteration
    
   // complex<double> omega = 1;
 
    //const double theta = (inverse ? -1 : 1) * 2 * M_PI / unityStep; // INVERSE
    complex<double> unityRoot = directionTwiddle[p-1] ;

    // each higher level doubles the step size
    for (int offset = 0; offset < N; offset += unityStep) {

   //    combine within a step segment (note only iterate over half step)
   danielson_loopunroll16(dual, directionTwiddle[p-1], unityStep/2, offset); 
  /*    for (int k = 0; k < unityStep/2; k++) {
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

void loopunrollFFT_8(const vector<complex<double> >& primal, vector<complex<double> >& dual, const vector<complex<double> >& directionTwiddle  ,const int P) 
{
  const int N = primal.size();
  const bool inverse = P < 0;
  
  const int absP = find_absP(P);	// Whether positive or negative, absP is always positive
  butterfly_loopunroll8(primal, dual, N, absP); 
  // bottom level of iteration tree --> puts elements in butterfly order
/* for (int i = 0; i < N; i++){ 
    dual[i] = primal[reverseBits(i) >> (32 - absP)];
  //  dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
   // dual[i+2] = primal[reverseBits(i+2) >> (32 - absP)];
    //dual[i+3] = primal[reverseBits(i+3) >> (32 - absP)];
  } 
 */ 
  
  //start of the danielson 
  // there are absP levels above the bottom
  for (int p = 1; p <= absP; p++) {

    // complex root of unity
    const int unityStep = 0x1 << p;	// --> starts at two, doubles each iteration
    
   // complex<double> omega = 1;
 
    //const double theta = (inverse ? -1 : 1) * 2 * M_PI / unityStep; // INVERSE
    complex<double> unityRoot = directionTwiddle[p-1] ;

    // each higher level doubles the step size
    for (int offset = 0; offset < N; offset += unityStep) {

   //    combine within a step segment (note only iterate over half step)
   danielson_loopunroll8(dual, directionTwiddle[p-1], unityStep/2, offset); 
  /*    for (int k = 0; k < unityStep/2; k++) {
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

void loopunrollFFT_4(const vector<complex<double> >& primal, vector<complex<double> >& dual, const vector<complex<double> >& directionTwiddle  ,const int P) 
{
  const int N = primal.size();
  const bool inverse = P < 0;
  
  const int absP = find_absP(P);	// Whether positive or negative, absP is always positive
  butterfly_loopunroll4(primal, dual, N, absP); 
  // bottom level of iteration tree --> puts elements in butterfly order
/* for (int i = 0; i < N; i++){ 
    dual[i] = primal[reverseBits(i) >> (32 - absP)];
  //  dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
   // dual[i+2] = primal[reverseBits(i+2) >> (32 - absP)];
    //dual[i+3] = primal[reverseBits(i+3) >> (32 - absP)];
  } 
 */ 
  
  //start of the danielson 
  // there are absP levels above the bottom
  for (int p = 1; p <= absP; p++) {

    // complex root of unity
    const int unityStep = 0x1 << p;	// --> starts at two, doubles each iteration
    
   // complex<double> omega = 1;
 
    //const double theta = (inverse ? -1 : 1) * 2 * M_PI / unityStep; // INVERSE
    complex<double> unityRoot = directionTwiddle[p-1] ;

    // each higher level doubles the step size
    for (int offset = 0; offset < N; offset += unityStep) {

   //    combine within a step segment (note only iterate over half step)
   danielson_loopunroll4(dual, directionTwiddle[p-1], unityStep/2, offset); 
  /*    for (int k = 0; k < unityStep/2; k++) {
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

void loopunrollFFT_2(const vector<complex<double> >& primal, vector<complex<double> >& dual, const vector<complex<double> >& directionTwiddle  ,const int P) 
{
  const int N = primal.size();
  const bool inverse = P < 0;
  
  const int absP = find_absP(P);	// Whether positive or negative, absP is always positive
  butterfly_loopunroll2(primal, dual, N, absP); 
  // bottom level of iteration tree --> puts elements in butterfly order
/* for (int i = 0; i < N; i++){ 
    dual[i] = primal[reverseBits(i) >> (32 - absP)];
  //  dual[i+1] = primal[reverseBits(i+1) >> (32 - absP)];
   // dual[i+2] = primal[reverseBits(i+2) >> (32 - absP)];
    //dual[i+3] = primal[reverseBits(i+3) >> (32 - absP)];
  } 
 */ 
  
  //start of the danielson 
  // there are absP levels above the bottom
  for (int p = 1; p <= absP; p++) {

    // complex root of unity
    const int unityStep = 0x1 << p;	// --> starts at two, doubles each iteration
    
   // complex<double> omega = 1;
 
    //const double theta = (inverse ? -1 : 1) * 2 * M_PI / unityStep; // INVERSE
    complex<double> unityRoot = directionTwiddle[p-1] ;

    // each higher level doubles the step size
    for (int offset = 0; offset < N; offset += unityStep) {

   //    combine within a step segment (note only iterate over half step)
   danielson_loopunroll2(dual, directionTwiddle[p-1], unityStep/2, offset); 
  /*    for (int k = 0; k < unityStep/2; k++) {
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
void basic_opt_FFT(const vector<complex<double> >& primal, vector<complex<double> >& dual,const int P) 
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
    for (int j = 0; j < N; j++)
      dual[j] /= N;
}

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


void timing(int n){
  struct timespec time1, time2, diffTime;
  struct timespec diff(struct timespec start, struct timespec end);  
  struct timespec diffArray[10]; 
  int clock_gettime(clockid_t clk_id, struct timespec *tp);

    



  for(int i = 1 ; i <= n; i++){ 
  	int N = 2; 
	N = pow(N,i);  
	int P = i;

  vector<complex<double> > primal(N, 0);
  for (int j = 0; j < N; j++)
    primal[j] = j;

  vector<complex<double> > dual0(N, 0);
  vector<complex<double> > dual1(N, 0);
  vector<complex<double> > dual2(N, 0);
  vector<complex<double> > dual3(N, 0);
  vector<complex<double> > dual4(N, 0);
  vector<complex<double> > dual5(N, 0);
  vector<complex<double> > dual6(N, 0);
  vector<complex<double> > dual7(N, 0);
  vector<complex<double> > dual8(N, 0);

  vector<complex<double> > forwardTwiddle;
  vector<complex<double> > backTwiddle;
  vector<complex<double> > dualPrime0(N, 0);
  vector<complex<double> > dualPrime1(N, 0);
  vector<complex<double> > dualPrime2(N, 0);
  vector<complex<double> > dualPrime3(N, 0);
  vector<complex<double> > dualPrime4(N, 0);
  vector<complex<double> > dualPrime5(N, 0);
  vector<complex<double> > dualPrime6(N, 0);
  vector<complex<double> > dualPrime7(N, 0);

  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
  twiddle_factors(forwardTwiddle, backTwiddle, N, P); 
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
  diffArray[0] = diff(time1, time2); 

  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
  baselineFFT(primal, dual0, P); // primal -> dual
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
  diffArray[1] = diff(time1, time2); 


  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
  basic_opt_FFT(primal, dual1, P); // primal -> dual
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
  diffArray[2] = diff(time1, time2); 


  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
  twiddleFFT(primal, dual2, forwardTwiddle, P); // primal -> dual
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
  diffArray[3] = diff(time1, time2); 


   
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
  loopunrollFFT_2(primal, dual3, forwardTwiddle, P); // primal -> dual
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
  diffArray[4] = diff(time1, time2); 
   
 
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
  loopunrollFFT_4(primal, dual4, forwardTwiddle, P); // primal -> dual
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
  diffArray[5] = diff(time1, time2); 
   
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
  loopunrollFFT_8(primal, dual5, forwardTwiddle, P); // primal -> dual
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
  diffArray[6] = diff(time1, time2); 
   
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
  loopunrollFFT_16(primal, dual6, forwardTwiddle, P); // primal -> dual
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
  diffArray[7] = diff(time1, time2); 
   
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
  loopunrollFFT_32(primal, dual7, forwardTwiddle, P); // primal -> dual
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
  diffArray[8] = diff(time1, time2); 
   
  // use -P as flag for inverse
   iterativeFFT(dual0, dualPrime0, -P); // dual -> primal
   iterativeFFT(dual1, dualPrime1, -P); // dual -> primal
   iterativeFFT(dual2, dualPrime2, -P); // dual -> primal 
   iterativeFFT(dual3, dualPrime3, -P); // dual -> primal
   iterativeFFT(dual4, dualPrime4, -P); // dual -> primal
   iterativeFFT(dual5, dualPrime5, -P); // dual -> primal
   iterativeFFT(dual6, dualPrime6, -P); // dual -> primal
   iterativeFFT(dual7, dualPrime7, -P); // dual -> primal
/*
    cout << N << "\t" << dualPrime0[N-1] << endl;
    cout << N << "\t" << dualPrime1[N-1] << endl;
    cout << N << "\t" << dualPrime2[N-1] << endl;
    cout << N << "\t" << dualPrime3[N-1] << endl;
  */
    dualPrime3.clear(); 



  printf("Numer of Elements:,  %d \n", N ); 
  printf("Twiddle Factor :,    %ld, cycles\n",(((long int)((double)(CPG)*(double)(GIG * (diffArray[0].tv_sec + diffArray[0].tv_nsec)))))); 
  printf("Baseline:,    %ld, cycles\n", (((long int)((double)(CPG)*(double)(GIG * (diffArray[1].tv_sec + diffArray[1].tv_nsec)))))); 
  printf("Basic Optimizations:,  %ld, cycles\n", (((long int)((double)(CPG)*(double)(GIG * (diffArray[2].tv_sec + diffArray[2].tv_nsec)))))); 
  printf("Twiddle Precomputations:,   %ld cycles\n", (((long int)((double)(CPG)*(double)(GIG * (diffArray[3].tv_sec + diffArray[3].tv_nsec)))))); 
  printf("Loop Unrolling 2:,  %ld, cycles\n", (((long int)((double)(CPG)*(double)(GIG * (diffArray[4].tv_sec + diffArray[4].tv_nsec)))))); 
  printf("Loop Unrolling 4:,  %ld, cycles\n", (((long int)((double)(CPG)*(double)(GIG * (diffArray[5].tv_sec + diffArray[5].tv_nsec)))))); 
  printf("Loop Unrolling 8:,  %ld, cycles\n", (((long int)((double)(CPG)*(double)(GIG * (diffArray[6].tv_sec + diffArray[6].tv_nsec)))))); 
  printf("Loop Unrolling 16:,  %ld, cycles\n", (((long int)((double)(CPG)*(double)(GIG * (diffArray[7].tv_sec + diffArray[7].tv_nsec)))))); 
  printf("Loop Unrolling 32:,  %ld, cycles\n", (((long int)((double)(CPG)*(double)(GIG * (diffArray[8].tv_sec + diffArray[8].tv_nsec)))))); 
  printf("\n"); 
  } 

}


int main(int argc, char *argv[]) {
  struct timespec time1, time2, diffTime;
  struct timespec diff(struct timespec start, struct timespec end);  
  struct timespec diffArray[5]; 
  int clock_gettime(clockid_t clk_id, struct timespec *tp);

  // input number of coefficients
  //cout << "input number of coefficients" << endl;
  
//  int N = 1048576*2*2*2*2; //use this for actual results 
 //int N = 16; 
  int N = 8; 
 
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

  vector<complex<double> > dual0(N, 0);
  vector<complex<double> > dual1(N, 0);
  vector<complex<double> > dual2(N, 0);
  vector<complex<double> > dual3(N, 0);

//timing the various optimizations 

  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
  (primal, dual0, P); // primal -> dual
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
  diffArray[0] = diff(time1, time2); 


  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
  basic_opt_FFT(primal, dual1, P); // primal -> dual
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
  diffArray[1] = diff(time1, time2); 


  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
  twiddleFFT(primal, dual2, forwardTwiddle, P); // primal -> dual
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
  diffArray[2] = diff(time1, time2); 


   
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
  loopunrollFFT_dep(primal, dual3, forwardTwiddle, P); // primal -> dual
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
  diffArray[3] = diff(time1, time2); 

 
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
  vector<complex<double> > dualPrime0(N, 0);
  vector<complex<double> > dualPrime1(N, 0);
  vector<complex<double> > dualPrime2(N, 0);
  vector<complex<double> > dualPrime3(N, 0);

  // use -P as flag for inverse
   iterativeFFT(dual0, dualPrime0, -P); // dual -> primal
   iterativeFFT(dual1, dualPrime1, -P); // dual -> primal
   iterativeFFT(dual2, dualPrime2, -P); // dual -> primal 
   iterativeFFT(dual3, dualPrime3, -P); // dual -> primal


 // cout << endl;
 //for (int i = 0; i < dualPrime3.size(); i++)
  //cout << i << "\t" << dualPrime3[i] << endl;

  // print dual coefficients
  cout << "dualPrime:" << endl;
 for (int i = 0; i < dual3.size(); i++){ 
    //cout << i << " dual0 " << "\t" << dual0[i] << " dual1 "  << "\t" << dual1[i] << " dual2 "<< "\t" << dual2[i] << " dual3 "<< "\t" << dual3[i] << endl;
  //  cout << "\t" << dual0[i] << "," <<   "\t"  << dual3[i] << endl;
   }
 
   // cout << N << "\t" << dualPrime0[N-1] << endl;
  //  cout << N << "\t" << dualPrime1[N-1] << endl;
 //   cout << N << "\t" << dualPrime2[N-1] << endl;
 //    cout << N << "\t" << dualPrime3[N-1] << endl;
  

//Timing output
//
  timing(24); 
/*
  diffTime = diff(time1, time2); 
  printf("Baseline:    %ld.%.9ld seconds\n", (long long)diffArray[0].tv_sec, diffArray[0].tv_nsec); 
  printf("Basic Optimizations:  %ld.%.9ld seconds\n", (long long)diffArray[1].tv_sec, diffArray[1].tv_nsec); 
  printf("Twiddle Precomputations:   %ld.%.9ld seconds\n", (long long)diffArray[2].tv_sec, diffArray[2].tv_nsec); 
  printf("222Loop Unrolling:  %ld.%.9ld seconds\n", (long long)diffArray[3].tv_sec, diffArray[3].tv_nsec); 
*/

  exit(0);
}
