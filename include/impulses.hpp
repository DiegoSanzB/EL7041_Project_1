#include <cmath>
#include <complex>
#include <valarray>
#include <unordered_map>
#include <functional>
#include <cassert>
#include <iostream>

#include "fft.hpp"

#ifndef IMPULSES
#define IMPULSES    

using namespace std;
const double lim = 1e-10;

typedef complex<double> Complex;
typedef valarray<Complex> CArray;

// PULSE SHAPES
double RC(double tau, double alpha);
inline function<double(double, unordered_map<string, double>)> RC_function = [](double tau, unordered_map<string, double> params) 
{
    return RC(tau, params["alpha"]);
};

double BTRC(double tau, double alpha);
inline function<double(double, unordered_map<string, double>)> BTRC_function = [](double tau, unordered_map<string, double> params)
{
    return BTRC(tau, params["alpha"]);
};

double IPLCP(double tau, double alpha, double mu, double epsilon, double gamma);
inline function<double(double, unordered_map<string, double>)> IPLCP_function = [](double tau, unordered_map<string, double> params)
{
    return IPLCP(tau, params["alpha"], params["mu"], params["epsilon"], params["gamma"]);
};

double ELP(double tau, double alpha, double betea);
inline function<double(double, unordered_map<string, double>)> ELP_function = [](double tau, unordered_map<string, double> params)
{
    return ELP(tau, params["alpha"], params["beta"]);
};

class Pulse {
    private:
        string pulse_type;
        function<double(double, unordered_map<string, double>)> pulse_function;
        unordered_map<string, double> pulse_params;
        // valarray<complex<double>> pulse; // TODO: maybe store the computed pulse
    public:
        Pulse(function<double(double, unordered_map<string, double>)> pulse,
            unordered_map<string, double> pulse_params, 
            string pulse_type)
            : pulse_function(pulse), pulse_params(pulse_params), pulse_type(pulse_type) {}
        string get_pulse_type();
        double evaluate(double, double);
        // double evaluate_truncated(double, double);
        valarray<double> get_array(int factor = 10, int fs = 10);
        valarray<double> get_array_fft(int factor = 10, int fs = 10, int NFFT = 1024);
};
#endif
// -- END OF FILE -- //
