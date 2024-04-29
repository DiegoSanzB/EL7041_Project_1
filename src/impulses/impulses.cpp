#include "impulses/impulses.hpp"


// void fft(CArray& x)
// {
//     const size_t N = x.size();
//     if (N <= 1) return;

//     // divide
//     CArray even = x[slice(0, N/2, 2)];
//     CArray  odd = x[slice(1, N/2, 2)];

//     // conquer
//     fft(even);
//     fft(odd);

//     // combine
//     for (size_t k = 0; k < N/2; ++k)
//     {
//         Complex t = polar(1.0, -2 * PI * k / N) * odd[k];
//         x[k    ] = even[k] + t;
//         x[k+N/2] = even[k] - t;
//     }
// }

// PULSES
// Raised Cosine
double RC(double tau, double alpha) {
    // Sinc
    double sincDen = PI * tau;
    double sinc = (abs(sincDen) > lim) ? sin(PI * tau) / sincDen : 1.0;
    // Cosine
    double cosDen = 1.0 - pow(2.0 * alpha * tau, 2.0);
    double cosine = (abs(cosDen) > lim) ? cos(alpha * PI * tau) / cosDen : 1.0;
    return sinc * cosine;
}

// Wrap it into a function RC_function
// function<double(double, unordered_map<string, double>)> RC_function = [](double tau, unordered_map<string, double> params) 
// {
//     return RC(tau, params["alpha"]);
// };


// Better than Raised Cosine Pulse (BTRC)
double BTRC(double tau, double alpha){
    // constants
    double beta = (PI * alpha)/log(2.0);
    // Sinc
    double sincDen = PI * tau;
    double sinc = (abs(sincDen) > lim) ? sin(PI * tau) / sincDen : 1.0;
    // The right term (does it have a name?)
    double right_Den = 1.0 + pow(beta * tau, 2.0);
    double right_Num = (((2.0 * beta * tau) * sin(PI * alpha * tau)) + (2.0 * cos(PI * alpha * tau)) - 1.0);
    double right = (abs(right_Den) > lim) ? right_Num / right_Den : 1.0;
    return sinc * right;
}

// Wrap it into a function BTRC_function
// inline function<double(double, unordered_map<string, double>)> BTRC_function = [](double tau, unordered_map<string, double> params) {
//         return BTRC(tau, params["alpha"]);
//     };


// Improved Parametric Linear Combination Pulse (IPLCP)
double IPLCP(double tau, double alpha, double mu, double epsilon, double gamma){
    // exp section
    double exp = std::exp((-1.0 * epsilon) * pow((PI * tau), 2.0));
    // Sinc
    double sincDen = PI * tau;
    double sinc = (abs(sincDen) > lim) ? sin(PI * tau) / sincDen : 1.0;
    // Section inside of a parenthesis
    double insideDen = pow((PI * alpha * tau), 2.0);
    double insideNum = (4.0 * (1.0 - mu) * sin(PI * alpha * tau * 0.5)) + (PI *  alpha * mu * tau * sin(PI * alpha * tau));
    double inside = (abs(insideDen) > lim) ? insideNum / insideDen : 1.0;
    double ans = exp * (pow((sinc * inside), gamma));
    return ans;
}


// Wrap it into a function IPLCP_function
// function<double(double, unordered_map<string, double>)> IPLCP_function = [](double tau, unordered_map<string, double> params) {
//         return IPLCP(tau, params["alpha"], params["mu"], params["epsilon"], params["gamma"]);
//     };


// Exponential Linear Pulse (ELP)
double ELP(double tau, double alpha, double beta){
    // exp section
    double exp = std::exp((-1.0 * PI) * (0.5 * beta) * pow(tau, 2.0));
    // Right sinc
    double RightSincDen = PI * tau;
    double RightSinc = (abs(RightSincDen) > lim) ? sin(PI * tau) / RightSincDen : 1.0;
    // Left sinc
    double LeftSincDen = PI * alpha * tau;
    double LeftSinc = (abs(LeftSincDen) > lim) ? sin(PI * alpha * tau) / LeftSincDen : 1.0;
    return exp * RightSinc * LeftSinc;
}


// Wrap it into a function ELP_function
// function<double(double, unordered_map<string, double>)> ELP_function = [](double tau, unordered_map<string, double> params) {
//         return ELP(tau, params["alpha"], params["beta"]);
//     };


// Shared class for pulses
// class Pulse {
// private:
//     string pulse_type;
//     function<double(double, unordered_map<string, double>)> pulse_function;
//     unordered_map<string, double> pulse_params;

// public:
    // Pulse class definition
    

// Pulse id
string Pulse::get_pulse_type()
{
    return pulse_type;
}

// Pulse evaluation in time domain
double Pulse::evaluate(double tau)
{
    return pulse_function(tau, pulse_params);
}

// Pulse array in time domain (say it's size in x*tau times, both sides)
valarray<double> Pulse::get_array(int factor, int fs)
{
    // Determine array size needed
    int arr_size = factor * 2 * fs + 1;
    // Create empty arrays (time and results holders)
    valarray<double> tau_array(0.0, arr_size);
    valarray<double> pulse_array(0.0, arr_size);
    // Fill time array
    for (int i = 0; i < arr_size; i++) {
        tau_array[i] = -factor + i / static_cast<double>(fs);   // TODO: check if the time domain is correct
    }
    // Now evaluate the pulse at each time
    for (int i = 0; i < arr_size; i++) {
        pulse_array[i] = evaluate(tau_array[i]);
    }
    return pulse_array;
}

// Pulse array in frequency domain
// #TODO: what if the input is the time domain pulse?
valarray<double> Pulse::get_array_fft(int factor, int fs, int NFFT)
{
    // Determine array size needed
    int arr_size = factor * 2 * fs + 1;
    // Assert NFFT > arr_size
    assert((NFFT > arr_size) && "NFFT must be greater than the pulse array size");
    // Create empty arrays (time and results holders)
    valarray<double> tau_array(0.0, arr_size);
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
    pulse_array_fft[slice(0, arr_size, 1)] = pulse_array;
    fft(pulse_array_fft);
    // Center the pulse
    pulse_array_fft = pulse_array_fft.cshift(NFFT / 2);
    // Abs and normalize
    pulse_array_fft = abs(pulse_array_fft)/static_cast<double>(fs);
    valarray<double> ans (NFFT);
    // Insert it into the valarray
    for (int i = 0; i < NFFT; i++) {
        ans[i] = pulse_array_fft[i].real();
    }
    return ans;
}
// };

// int main() {
//     // Define the pulse
//     Pulse pulse(RC_function, {{"alpha", 0.35}}, "RC");
//     // CArray pulse_array = pulse.get_array();
//     // // print the pulse array
//     // for (int i = 0; i < pulse_array.size(); i++) {
//     //     cout << pulse_array[i].real() << ',';
//     // }
//     // Print the fft
//     valarray<double> pulse_array_fft = pulse.get_array_fft();
//     for (int i = 0; i < pulse_array_fft.size(); i++) {
//         cout << pulse_array_fft[i] << ',';
//     }

//     cout << endl;
//     return 0;
// }

