#include <cmath>
#include <valarray>
#include <random>

#include "fft.hpp"
#include "impulses/impulses.hpp"

const int time_unit_period = 1;
const valarray<double> OFFSETS = {0.05, 0.1, 0.2, 0.25};

valarray<double> ber_isi(Pulse pulse, double fs, double snr, int nbits, int M);
valarray<double> ber_cci(Pulse pulse, double fs, double snr, double sir, int L, int nbits, int M);
valarray<double> ber_isi_cci(Pulse pulse, double fs, double snr, double sir, int L, int nbits, int M);

// -- END OF FILE -- //