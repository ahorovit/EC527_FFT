#pragma once

#include <pthread.h>
#include <unistd.h>
#include <queue>
#include <cstdlib>

// For print statements, paste:
//		#define THREADPOOL_ENABLE_PRINT
// before:
//		#include "thread_pool.hpp"
// In files that include this one.
#ifdef THREADPOOL_ENABLE_PRINT
	#include <iostream>
#endif

void * aquireWork(void * threadarg);

class ThreadPool {
private:
	// Vector of threads stored in the pool
	std::vector<pthread_t> pool;
		
	// Begins thread execution
	void startThreads();	
	
public:
	// Default constructor, populates pool with number of threads equal to number of CPU cores
	ThreadPool();
	
	// Destructor, rejoins threads when an instance of ThreadPool goes out of scope
	~ThreadPool();
	
	// Returns true if there is work remaining
	bool hasWork();
	
	// Populates pool with number of threads equal to the value of the input parameter
	ThreadPool(std::size_t);
	
	// Work functions allowed are of the same format as those used by pthreads
	void addWork(void * (*)(void *), void *);
};
