// #include "impulses/impulses.hpp"
#include <cmath>
#include <complex>
#include <valarray>
#include <unordered_map>
#include <functional>
#include <cassert>
#include <iostream>

const double PI = 3.141592653589793238460;
const double lim = 1e-10;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

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

// PULSES
// Raised Cosine
double RC(double tau, double alpha) {
    // Sinc
    double sincDen = PI * tau;
    double sinc = (std::abs(sincDen) > lim) ? std::sin(PI * tau) / sincDen : 1.0;
    // Cosine
    double cosDen = 1.0 - std::pow(2.0 * alpha * tau, 2.0);
    double cosine = (std::abs(cosDen) > lim) ? std::cos(alpha * PI * tau) / cosDen : 1.0;
    return sinc * cosine;
}
// Wrap it into a function RC_function
std::function<double(double, std::unordered_map<std::string, double>)> RC_function = [](double tau, std::unordered_map<std::string, double> params) {
        return RC(tau, params["alpha"]);
    };

// Better than Raised Cosine Pulse (BTRC)
double BTRC(double tau, double alpha){
    // constants
    double beta = (PI * alpha)/std::log(2.0);
    // Sinc
    double sincDen = PI * tau;
    double sinc = (std::abs(sincDen) > lim) ? std::sin(PI * tau) / sincDen : 1.0;
    // The right term (does it have a name?)
    double right_Den = 1.0 + std::pow(beta * tau, 2.0);
    double right_Num = (((2.0 * beta * tau) * std::sin(PI * alpha * tau)) + (2.0 * std::cos(PI * alpha * tau)) - 1.0);
    double right = (std::abs(right_Den) > lim) ? right_Num / right_Den : 1.0;
    return sinc * right;
}
// Wrap it into a function BTRC_function
std::function<double(double, std::unordered_map<std::string, double>)> BTRC_function = [](double tau, std::unordered_map<std::string, double> params) {
        return BTRC(tau, params["alpha"]);
    };

// Improved Parametric Linear Combination Pulse (IPLCP)
double IPLCP(double tau, double alpha, double mu, double epsilon, double gamma){
    // exp section
    double exp = std::exp((-1.0 * epsilon) * std::pow((PI * tau), 2.0));
    // Sinc
    double sincDen = PI * tau;
    double sinc = (std::abs(sincDen) > lim) ? std::sin(PI * tau) / sincDen : 1.0;
    // Section inside of a parenthesis
    double insideDen = std::pow((PI * alpha * tau), 2.0);
    double insideNum = (4.0 * (1.0 - mu) * std::sin(PI * alpha * tau * 0.5)) + (PI *  alpha * mu * tau * std::sin(PI * alpha * tau));
    double inside = (std::abs(insideDen) > lim) ? insideNum / insideDen : 1.0;
    double ans = exp * (std::pow((sinc * inside), gamma));
    return ans;
}
// Wrap it into a function IPLCP_function
std::function<double(double, std::unordered_map<std::string, double>)> IPLCP_function = [](double tau, std::unordered_map<std::string, double> params) {
        return IPLCP(tau, params["alpha"], params["mu"], params["epsilon"], params["gamma"]);
    };

// Exponential Linear Pulse (ELP)
double ELP(double tau, double alpha, double beta){
    // exp section
    double exp = std::exp((-1.0 * PI) * (0.5 * beta) * std::pow(tau, 2.0));
    // Right sinc
    double RightSincDen = PI * tau;
    double RightSinc = (std::abs(RightSincDen) > lim) ? std::sin(PI * tau) / RightSincDen : 1.0;
    // Left sinc
    double LeftSincDen = PI * alpha * tau;
    double LeftSinc = (std::abs(LeftSincDen) > lim) ? std::sin(PI * alpha * tau) / LeftSincDen : 1.0;
    return exp * RightSinc * LeftSinc;
}
// Wrap it into a function ELP_function
std::function<double(double, std::unordered_map<std::string, double>)> ELP_function = [](double tau, std::unordered_map<std::string, double> params) {
        return ELP(tau, params["alpha"], params["beta"]);
    };

// Shared class for pulses
class Pulse {
private:
    std::string pulse_type;
    std::function<double(double, std::unordered_map<std::string, double>)> pulse_function;
    std::unordered_map<std::string, double> pulse_params;

public:
    // Pulse class definition
    Pulse(std::function<double(double, std::unordered_map<std::string, double>)> pulse,
          std::unordered_map<std::string, double> pulse_params, std::string pulse_type)
          : pulse_function(pulse), pulse_params(pulse_params), pulse_type(pulse_type) {}

    // Pulse id
    std::string get_pulse_type() const {
        return pulse_type;
    }

    // Pulse evaluation in time domain
    double evaluate(double tau) const {
        return pulse_function(tau, pulse_params);
    }

    // Pulse array in time domain (say it's size in x*tau times, both sides)
    std::valarray<double> get_array(int factor = 10, int fs = 10) const {
        // Determine array size needed
        int arr_size = factor * 2 * fs + 1;
        // Create empty arrays (time and results holders)
        std::valarray<double> tau_array(0.0, arr_size);
        std::valarray<double> pulse_array(0.0, arr_size);
        // Fill time array
        for (int i = 0; i < arr_size; i++) {
            tau_array[i] = -fs + i / static_cast<double>(fs);
        }
        // Now evaluate the pulse at each time
        for (int i = 0; i < arr_size; i++) {
            pulse_array[i] = evaluate(tau_array[i]);
        }
        return pulse_array;
    }

    // Pulse array in frequency domain
    // #TODO: what if the input is the time domain pulse?
    std::valarray<double> get_array_fft(int factor = 10, int fs = 10, int NFFT = 1024) const {
        // Determine array size needed
        int arr_size = factor * 2 * fs + 1;
        // Assert NFFT > arr_size
        assert((NFFT > arr_size) && "NFFT must be greater than the pulse array size");
        // Create empty arrays (time and results holders)
        std::valarray<double> tau_array(0.0, arr_size);
        CArray pulse_array(0.0, arr_size);
        // Fill time array
        for (int i = 0; i < arr_size; i++) {
            tau_array[i] = -fs + i / static_cast<double>(fs);
        }
        // Now evaluate the pulse at each time
        for (int i = 0; i < arr_size; i++) {
            pulse_array[i] = evaluate(tau_array[i]);
        }
        // Now do the FFT
        CArray pulse_array_fft(NFFT);
        pulse_array_fft = 0.0; // Fill with zeros
        pulse_array_fft[std::slice(0, arr_size, 1)] = pulse_array;
        fft(pulse_array_fft);
        // Center the pulse
        pulse_array_fft = pulse_array_fft.cshift(NFFT / 2);
        // Abs and normalize
        pulse_array_fft = std::abs(pulse_array_fft)/static_cast<double>(fs);
        std::valarray<double> ans (NFFT);
        // Insert it into the valarray
        for (int i = 0; i < NFFT; i++) {
            ans[i] = pulse_array_fft[i].real();
        }
        return ans;
    }
};

int main() {
    // Define the pulse
    Pulse pulse(RC_function, {{"alpha", 0.35}}, "RC");
    // CArray pulse_array = pulse.get_array();
    // // print the pulse array
    // for (int i = 0; i < pulse_array.size(); i++) {
    //     std::cout << pulse_array[i].real() << ',';
    // }
    // Print the fft
    std::valarray<double> pulse_array_fft = pulse.get_array_fft();
    for (int i = 0; i < pulse_array_fft.size(); i++) {
        std::cout << pulse_array_fft[i] << ',';
    }

    std::cout << std::endl;
    return 0;
}