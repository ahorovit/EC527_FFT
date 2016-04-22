#include "fft.hpp"

namespace fft527 {
	/*void fastFourierTransform(complexVector & elements, float dir) {
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
	}*/
	
	void fastFourierTransform(complexVector & elements) {
		size_t length = elements.size();
		unsigned int k, n;
		k = length;
		double thetaT = PI / length;
		complex phiT = complex(cos(thetaT), sin(thetaT));
		complex T;
		
		while (k > 1) {
			n = k;
			k >>= 1;
			phiT = phiT * phiT;
			T = 1.0L;
			for (unsigned int l = 0; l < k; l++) {
				for (unsigned int a = l; a < length; a += n) {
					unsigned int b = a + k;
					complex t = elements[a] - elements[b];
					elements[a] += elements[b];
					elements[b] = t * T;
				}
				T *= phiT;
			}
		}
		
		unsigned int m = (unsigned int) log2(length);
		for (unsigned int a = 0; a < length; a++) {
			unsigned int b = a;
			b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
			b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
			b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
			b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
			b = ((b >> 16) | (b << 16)) >> (32 - m);
			if (b > a) {
				complex t = elements[a];
				elements[a] = elements[b];
				elements[b] = t;
			}
		}
	}

	void fftForward(complexVector & elements) {
		fastFourierTransform(elements);
	}

	void fftReverse(complexVector & elements) {
		for (int i = 0; i < elements.size(); i++) {
			elements[i] = std::conj(elements[i]);
		}
		
		fastFourierTransform(elements);
		
		size_t length = elements.size();
		for (int i = 0; i < length; i++) {
			elements[i] = std::conj(elements[i]);
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