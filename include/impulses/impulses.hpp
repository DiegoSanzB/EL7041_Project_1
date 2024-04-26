#include <complex>
#include <valarray>

#include "fft.hpp"

const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

// PULSES
CArray raised_cosine(float alpha);
CArray raised_cosine(float alpha, float fs, int NFFT);

// CArray BTRC(float alpha, float beta);

// CArray IPLCP(float alpha, float mu, float epsilon, float gamma);
typedef std::valarray<Complex> CArray;
// CArray ELP(float alpha, float beta);

class Pulses {
    private:
        string pulse_type;
        valarray<complex<double>> pulse;
    public:

        
};

