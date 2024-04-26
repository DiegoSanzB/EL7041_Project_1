#include <valarray>

#include "fft.hpp"
#include "impulses/impulses.hpp"

const int unit_time = 1;
const valarray<float> OFFSETS = {0.05, 0.1, 0.2, 0.25};

CArray ber(CArray &pulse, int snr, int nbits, int M, float omega);

// -- END OF FILE -- //