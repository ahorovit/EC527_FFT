#include "fft.hpp"

namespace fft527 {
	void fastFourierTransform(complexVector & elements, float dir) {
		size_t length = elements.size();
		if (length <= 1)
			return;
			
		// Divide into even and odd halves
		complexVector evenHalf, oddHalf;
		evenHalf.reserve(length / 2);
		oddHalf.reserve(length / 2);
		for (int i = 0; i < length; i++) {
			if (i % 2 == 0) {
				evenHalf.push_back(elements[i]);
			} else { 
				oddHalf.push_back(elements[i]);
			}
		}
		
		// Conquer
		fastFourierTransform(evenHalf, dir);
		fastFourierTransform(oddHalf, dir);
		
		for (int i = 0; i < length / 2; ++i) {
			complex t = std::polar(1.0, (dir * 2) * PI * i / length) * oddHalf[i];
			elements[i] = evenHalf[i] + t;
			elements[i + length / 2] = evenHalf[i] - t;
		}
	}

	void fftForward(complexVector & elements) {
		fastFourierTransform(elements, -1.f);
	}

	void fftReverse(complexVector & elements) {
		fastFourierTransform(elements, 1.f);
		// Rescale
		size_t length = elements.size();
		for (int i = 0; i < length; i++) {
			elements[i] /= length;
		}
	}
	
	void dispReal(complexVector & elements) {
		for (auto & element : elements) {
			std::cout << element.real() << std::endl;
		}
	}
	
	void dispComplex(complexVector & elements) {
		for (auto & element : elements) {
			std::cout << element << std::endl;
		}
	}
}