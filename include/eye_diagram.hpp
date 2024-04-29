#include <random>

#include "fft.hpp"
#include "impulses/impulses.hpp"

// TODO: change pulse to Pulse class maybe
valarray<double> eye_diagram(size_t N, valarray<double> pulse, float fs);

valarray<int> random_bpsk(size_t N);

valarray<int> bpsk_upsample(valarray<int> &bpsk, size_t N, float fs);

valarray<double> conv(valarray<double> &a, valarray<double> &b);

// -- END OF FILE -- //