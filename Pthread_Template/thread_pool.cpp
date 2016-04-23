#include "thread_pool.hpp"

std::queue<std::pair<void * (*)(void *), void *>> work;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
bool terminate = false;
// Count of number of tasks enqueued, incremented when adding work and decremented upon completion
long int workCount = 0;

void * aquireWork(void * tid) {
	bool hasWork = false;
	void * (*workFn)(void *);
	void * workData;
	// Loop runs until user calls ThreadPool::releasePool(), or the process ends
	while(true) {
		if (terminate) {
			pthread_exit(nullptr);
		}
		
		// Wait for access to lock before interacting with the queue
		pthread_mutex_lock(&mutex);
		if (!work.empty()) {
			// Take work off the queue
			workFn = work.front().first;
			workData = work.front().second;
			work.pop();
			hasWork = true;
		}
		pthread_mutex_unlock(&mutex);
		
		// Execute work removed from the queue
		if (hasWork) {
			workFn(workData);
			--workCount;
			hasWork = false;
		}
	}
}
	
ThreadPool::ThreadPool() {
	pthread_t thread;
	int numCores = sysconf(_SC_NPROCESSORS_ONLN);
	pool.reserve(numCores);
	for (int i = 0; i < numCores; i++) {
		pool.push_back(thread);
	}
	startThreads();
}

ThreadPool::ThreadPool(std::size_t numThreads) {
	pthread_t thread;
	pool.reserve(numThreads);
	for (int i = 0; i < numThreads; i++) {
		pool.push_back(thread);
	}
	startThreads();
}

void ThreadPool::startThreads() {
	int rc, tid = 0;
	for (auto & element : pool) {
		rc = pthread_create(&element, nullptr, aquireWork, nullptr);
		if (rc) {
			#ifdef THREADPOOL_ENABLE_PRINT
				std::cerr << "pthread_create() failed with exit code: " << rc << std::endl;
			#endif
			exit(EXIT_FAILURE);
		}
	}
}

void ThreadPool::addWork(void * (* workFn)(void *), void * workData) {
	pthread_mutex_lock(&mutex);
	work.push(std::make_pair(workFn, workData));
	pthread_mutex_unlock(&mutex);
	workCount++;
}

ThreadPool::~ThreadPool() {
	terminate = true;
	size_t poolSize = pool.size();
	int rc;
	for (int i = 0; i < poolSize; i++) {
		rc = pthread_join(pool[i], nullptr);
		if (rc) {
			#ifdef THREADPOOL_ENABLE_PRINT
				std::cerr << "pthread_join() failed with exit code: " << rc << std::endl;
			#endif
			exit(EXIT_FAILURE);
		}
	}
	#ifdef THREADPOOL_ENABLE_PRINT
		std::cout << "Terminate succeded!" << std::endl;
	#endif
}

bool ThreadPool::hasWork() {
	return workCount > 0;
}
