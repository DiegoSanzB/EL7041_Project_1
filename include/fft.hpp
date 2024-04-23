#include <complex>
#include <valarray>

#ifndef FFT
#define FFT

using namespace std;

const double PI = 3.141592653589793238460;

typedef complex<double> Complex;
typedef valarray<Complex> CArray;
typedef valarray<double> RArray;

void fft_inplace(CArray&);
CArray fft(CArray&);

#endif

// -- END OF FILE -- // 