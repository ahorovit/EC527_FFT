all: testThreadPool.o thread_pool.o
	clang++ -std=c++11 -O2 testThreadPool.o thread_pool.o -o testPool.out

testThreadPool.o:
	clang++ -std=c++11 -O2 -c testThreadPool.cpp

thread_pool.o:
	clang++ -std=c++11 -O2 -c thread_pool.cpp

clean:
	rm testThreadPool.o thread_pool.o

cleanest:
	rm testThreadPool.o thread_pool.o testPool.out

