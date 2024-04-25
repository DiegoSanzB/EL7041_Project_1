#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
#include <valarray>

const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

// Define constant variables
const double fs = 10; // Sampling frequency
const int NFFT = 1024; // Number of points in the FFT
const double lim = 1e-10; // Limit for the denominator 

CArray BTRC(double alpha, bool time = true){
    // Generate Better Than Raised Cosine pulse
    // alpha: roll-off factor
    // time: if true, return the time domain pulse
    CArray FPulse (NFFT);
    for (int i = 0; i < NFFT; i++) {
        FPulse[i] = -0.5 + (i)/(NFFT);
    }
    CArray SampledPulse (2*fs*fs + 1);
    for (int i = 0; i <= 2*fs*fs; i++) {
        SampledPulse[i] = -fs + i/fs;
    }
    float beta = PI * alpha / log(2.0);
    // Sinc section
    CArray sincNum = sin(PI*SampledPulse);
    CArray sincDen = PI*SampledPulse;
    CArray sincOp (2*fs*fs + 1);
    // Replace with 1 if the value is less than the limit
    for (int i = 0; i <= 2*fs*fs; i++) {
        if (std::abs(sincDen[i]) > lim)
        {
            sincOp[i] = sincNum[i] / sincDen[i];
        }
        else {
            sincOp[i] = 1;
        }
    }
    // The other section
}