#pragma once

template <typename T, std::size_t S> 
struct complex_array {
	// Data field
	std::complex<T> data[S ? S : 1];
	
	// Iterators
	std::complex<T> * begin() {
		return addressof(data[0]);
	}
	
	std::complex<T> * end() {
		return addressof(data[S - 1]);
	}
	
	// Const Iterators
	std::complex<T> * begin() const {
		return addressof(data[0]);
	}
	
	std::complex<T> * end() const {
		return addressof(data[S - 1]);
	}
	
	// Subscript operator
	std::complex<T> & operator[](std::size_t idx) {
		return data[idx];
	}
	
	std::complex<T> & operator[](std::size_t idx) const {
		return data[idx];
	}
	
	// Size accessor
	constexpr std::size_t size() const {
		return S;
	}
};
