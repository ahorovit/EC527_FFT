all: test.o fft.o
	clang++ -std=c++11 test.o fft.o

test.o: fft.o
	clang++ -std=c++11 -c test.cpp

fft.o:
	clang++ -std=c++11 -c fft.cpp

clean:
	rm fft.o test.o
