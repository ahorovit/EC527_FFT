#include <iostream>
#include "fft.hpp"
#define SIZE 32

int main(int argc, char *argv[]) {
	fft527::complex test[SIZE];
	for (int i = 0; i < SIZE; i++) {
		test[i] = 1.0 * i;
	}
	fft527::complexVector testVec(test, test + SIZE);
	
	fft527::ComplexArray<int, 8> arr;
	
	fft527::dispReal(testVec);
	
	fft527::fftForward(testVec);
	fft527::dispComplex(testVec);
	
	fft527::fftReverse(testVec);
	fft527::dispReal(testVec);
	
	for (auto & element : arr) {
		
	}
}
