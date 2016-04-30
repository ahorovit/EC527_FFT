#include <climits>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <time.h> 

uint32_t reverseBits(uint32_t i) {
  uint32_t mask = 0x55555555; // 0101...
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
inline void iterativeFFT(const std::vector<std::complex<double> >& primal, std::vector<std::complex<double> >& dual, const int P) {
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
 
    const double theta = (is_inverse(inverse)) * 2 * M_PI / unityStep; // INVERSE
    const std::complex<double> unityRoot(cos(theta), sin(theta));

    // each higher level doubles the step size
    for (int offset = 0; offset < primal.size(); offset += unityStep) {
      std::complex<double> omega = 1;

      // combine within a step segment (note only iterate over half step)
      for (int k = 0; k < unityStep/2; k++) {
        const std::complex<double> u = dual[offset + k];

        const std::complex<double> t = omega * dual[offset + k + unityStep/2];
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