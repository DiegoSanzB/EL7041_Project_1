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

// Cooley–Tukey FFT (in-place, divide-and-conquer)
// Higher memory requirements and redundancy although more intuitive
void fft(CArray& x)
{
    const size_t N = x.size();
    if (N <= 1) return;

    // divide
    CArray even = x[std::slice(0, N/2, 2)];
    CArray  odd = x[std::slice(1, N/2, 2)];

    // conquer
    fft(even);
    fft(odd);

    // combine
    for (size_t k = 0; k < N/2; ++k)
    {
        Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        x[k    ] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }
}

CArray RC(double alpha, bool time = true){
    // Generate the raised cosine pulse
    // alpha: roll-off factor
    // time: if true, return the time domain pulse
    //       if false, return the frequency domain pulse
    CArray FPulse (NFFT);
    for (int i = 0; i < NFFT; i++) {
        FPulse[i] = -0.5 + (i)/(NFFT);
    }
    CArray SampledPulse (2*fs*fs + 1);
    for (int i = 0; i <= 2*fs*fs; i++) {
        SampledPulse[i] = -fs + i/fs;
    }
    // Sinc section of the pulse
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
    // Cosine section of the pulse
    CArray cosNum = cos(alpha*PI*SampledPulse);
    CArray cosDen = 1 - (2*alpha*SampledPulse)*(2*alpha*SampledPulse);
    CArray cosOp (2*fs*fs + 1);
    // Replace with 1 if the value is less than the limit
    for (int i = 0; i <= 2*fs*fs; i++) {
        if (std::abs(cosDen[i]) > lim)
        {
            cosOp[i] = cosNum[i] / cosDen[i];
        }
        else {
            cosOp[i] = 1;
        }
    }
    // Raised cosine pulse
    // Multiply the sinc and cosine sections
    CArray RCArray (2*fs*fs + 1);
    for (int i = 0; i <= 2*fs*fs; i++) {
        RCArray[i] = sincOp[i]*cosOp[i];
        // std::cout << RCArray[i].real() << ", ";
    }
    if (time) { 
        // Return the time domain pulse
        return RCArray;
    }
    else {
        // Return the frequency domain pulse
        CArray RCArrayFFT (NFFT);
        RCArrayFFT[std::slice(0, (2*fs*fs + 1), 1)] = RCArray;
        fft(RCArrayFFT);
        // Center the pulse
        RCArrayFFT = RCArrayFFT.cshift((NFFT + 2*fs*fs + 1)/2);
        return RCArrayFFT;
    }
    
}

int main() {
    double alpha = 0.35;
    // std::cout << "Raised cosine pulse with roll-off factor: " << alpha << std::endl;
    CArray RCArray = RC(alpha, false);
    for (int i = 0; i < RCArray.size(); i++){
        std::cout << std::abs(RCArray[i]) << ", ";
    }
    return 0;
}