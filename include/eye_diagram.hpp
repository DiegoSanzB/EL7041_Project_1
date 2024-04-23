#include <random>

#include "fft.hpp"

CArray eye_diagram(size_t N, CArray &pulse, float fs);

valarray<int> random_bpsk(size_t N);

valarray<int> bpsk_upsample(valarray<int> &bpsk, size_t N, float fs);

CArray conv(CArray &a, CArray &b);

// -- END OF FILE -- //