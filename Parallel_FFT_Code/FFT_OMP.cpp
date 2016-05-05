#include "FFT_equil.hpp"
#include <omp.h>
#include <vector>
#include <array>
#include <complex>
#include <iostream>
#include <time.h>

#define MAT_WIDTH 8192

// To Compile: g++ -Wall -fopenmp -std=c++0x FFT_OMP.cpp -o FFT_OMP
// For timing: time ./FFT_OMP

int main(void) {
	omp_set_num_threads(4);
	std::array<std::vector<std::complex<double> >, MAT_WIDTH> primals, duals;//, primalBuffer;

	const int P = lg(MAT_WIDTH);

	for (int i = 0; i < MAT_WIDTH; i++) {
		primals[i].resize(MAT_WIDTH);
		duals[i].resize(MAT_WIDTH);
		//primalBuffer[i].resize(MAT_WIDTH);
		for (int j = 0; j < MAT_WIDTH; j++) {
			primals[i][j] = j;
		}
	}

	#pragma omp parallel for
	for (int i = 0; i < MAT_WIDTH; i++) {
		iterativeFFT(primals[i], duals[i], P);
	}

	// Rotate vectors, for less false sharing within the column-wise FFTs. Begin with a transpose...
   	for (int i = 0; i < MAT_WIDTH; i++) {
        	for (int j = i + 1; j < MAT_WIDTH; j++) {
			primals[j][i] = duals[i][j];
        	}
    	}

	// then swap the columns.
	std::complex<double> tmp;
    	for (int i = 0; i < MAT_WIDTH; i++) {
        	for (int j = 0; j < MAT_WIDTH / 2; j++) {
            		tmp = primals[i][j];
            		primals[i][j] = primals[i][MAT_WIDTH - 1 - j];
            		primals[i][MAT_WIDTH - 1 - j] = tmp;
        	}
    	}

	#pragma omp parallel for
	for (int i = 0; i < MAT_WIDTH; i++) {
		iterativeFFT(primals[i], duals[i], P);
	}

	return EXIT_SUCCESS;
}

