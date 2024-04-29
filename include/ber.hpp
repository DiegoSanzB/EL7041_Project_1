#include <valarray>
#include <random>

#include "fft.hpp"
#include "impulses/impulses.hpp"

const int time_unit_period = 1;
const valarray<double> OFFSETS = {0.05, 0.1, 0.2, 0.25};

valarray<double> ber_isi(CArray &pulse, int snr, int nbits, int M, float omega);

// -- END OF FILE -- //