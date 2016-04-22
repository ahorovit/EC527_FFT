#include <vector>
#include <complex>
#include <iostream>

namespace fft527 {
	//#include "complex_array.h"
	#define PI 3.141592653589793238460
	#define data_t double

	typedef std::complex<data_t> complex;
	typedef std::vector<complex> complexVector;
	typedef std::vector<std::vector<complex>> complexMat;

	// A bi-directional single dimensional fft, and a bi-directional parameterized multidimensional fft
	// The wrapper functions below call these, you should not need to call them explicitly
	void fastFourierTransform(complexVector &, float);


	// A wrapper function that correctly calls fastFourierTransform(...) in the forward direction
	// Arguments:
	// 	arg0: Lvalue reference to a complexVector
	void fftForward(complexVector &);


	// Another wripper function that correctly calls fastFourierTransform(...) in the reverse direction
	// Arguments:
	//  arg0: Lvalue reference to a complexVector
	void fftReverse(complexVector &);

	void showReal(complexVector &);
	
	void dispReal(complexVector &);
	
	void dispReal(complexMat &);
	
	void dispComplex(complexVector &);
	
	void dispComplex(complexMat &);
}
