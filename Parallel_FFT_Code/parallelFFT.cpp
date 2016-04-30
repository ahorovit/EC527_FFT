#include <pthread.h>
#include <array>
#include <iostream>
#include <vector>
#include <complex>
#include "FFT_equil.hpp"

#define NUM_THREADS 2
#define VEC_LENGTH 32768

struct ThreadData {
	std::vector<std::complex<double> > * pPrimal;
	std::vector<std::complex<double> > * pDual;
	int P;
};

void * fftHelper(void * threadArg) {
	ThreadData td = *(static_cast<ThreadData *> (threadArg));
	iterativeFFT(*td.pPrimal, *td.pDual, td.P);
	pthread_exit(nullptr);
}

int main(int argc, char *argv[]) {
	std::array<pthread_t, NUM_THREADS> threads;
	std::array<ThreadData, NUM_THREADS> threadData;
	
	const int P = lg(VEC_LENGTH);
	
	std::array<std::vector<std::complex<double> >, NUM_THREADS> primals;
	std::array<std::vector<std::complex<double> >, NUM_THREADS> duals;
	std::array<std::vector<std::complex<double> >, NUM_THREADS> dualPrimes;
	for (int i = 0; i < NUM_THREADS; i++) {
		primals[i].resize(VEC_LENGTH);
		dualPrimes[i].resize(VEC_LENGTH);
		duals[i].resize(VEC_LENGTH);
		for (int j = 0; j < VEC_LENGTH; j++) {
			primals[i][j] = j;
		}
	}
	
	for (int i = 0; i < NUM_THREADS; i++) {
		threadData[i].pPrimal = &primals[i];
		threadData[i].pDual = &duals[i];
		threadData[i].P = P;
		int rc = pthread_create(&threads[i], nullptr, fftHelper, static_cast<void *>(&threadData[i]));
		if (rc)
			std::cerr << "pthread_create failed with error code: " << rc << std::endl;
	}
		
	for (auto it = threads.begin(); it != threads.end(); ++it) {
		int rc = pthread_join(*it, nullptr);
		if (rc)
			std::cerr << "pthread_join failed with error code: " << rc << std::endl;
	}
	
	std::cout << duals[0][VEC_LENGTH - 1] << " " << duals[0][VEC_LENGTH - 1] << std::endl;
	
	return EXIT_SUCCESS;
}
