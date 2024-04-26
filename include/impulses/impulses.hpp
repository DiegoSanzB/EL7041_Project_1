#include <complex>
#include <valarray>

#include "fft.hpp"

// PULSES
CArray raised_cosine(float alpha);
CArray raised_cosine(float alpha, float fs, int NFFT);
float raised_cosine_t(float unit_time, float t, float alpha);

CArray BTRC(float alpha, float beta);

CArray IPLCP(float alpha, float mu, float epsilon, float gamma);

CArray ELP(float alpha, float beta);

