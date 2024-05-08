#include <cmath>
#include <valarray>
#include <random>

#include "impulses.hpp"

const int time_unit_period = 1;
const valarray<double> OFFSETS = {0.05, 0.1, 0.2, 0.25};

valarray<double> ber_isi(Pulse pulse, double omega, double snr, int nbits, int M, double truncation);
valarray<double> ber_cci(Pulse pulse, double omega, double snr, double sir, int L, int nbits, int M, double truncation);
valarray<double> ber_isi_cci(Pulse pulse, double omega, double snr, double sir, int L, int nbits, int M, double truncation);

// -- END OF FILE -- //