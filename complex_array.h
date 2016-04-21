#pragma once
#include <iterator>

template <typename T, std::size_t S> 
struct ComplexArray {
	std::complex<T> data[S];
	std::size_t size = S;
	
	std::complex<T> * begin() { return addressof(data[0]); }
	
	std::complex<T> * end() { return addressof(data[S]); }
	
};
