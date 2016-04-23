//
//	Test file for demonstrating that the ThreadPool works (and simplifies working with pthreads).
//	thread_pool.hpp defines the class ThreadPool, which has the following functionality:
//		The default constructor creates a pool containing a number of threads equal to the number
//		of detected processor cores. You can optionally construct an instance of a ThreadPool with
//		a parameter for the number of desired threads.
//		
//		In order to add work to an instance of ThreadPool, simply call the addWork() function.
//		addWork takes as its first argument a function of type 'void * (*)(void *)' and a parameter
//		to pass to that function of type 'void *'.
//

#include "thread_pool.hpp"
#include <iostream>

int testValue = 0;

void * incrementFn(void * data) {
	testValue++;
}

struct workdata {
	// Whatever data you'd like to be passed to your function
};

int main(int argc, char *argv[]) {
	ThreadPool testPool;
	
	// Create a thread data argument for the function
	workdata wd;
	void * workData = static_cast<void *> (&wd);
	
	// Add some test work
	for (int i = 0; i < 10; i++)
	testPool.addWork(incrementFn, workData);

	// Wait for it to finish its work
	while (testPool.hasWork());
	
	// Demonstrate that it works
	std::cout << testValue << std::endl;
	
	return EXIT_SUCCESS;
}
