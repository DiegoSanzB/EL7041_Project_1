#include <cmath>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <string>

#include "ber.hpp"
#include "data_manager.hpp"
#include "eye_diagram.hpp"
#include "fft.hpp"
#include "impulses/impulses.hpp"

using namespace std;

int main() {  
    double alpha = 0.35, omega = 0.10;
    int factor = 10, fs = 10, NFFT = pow(2, fs);
    int snr_db = 15, sir_db = 10, L = 2, nbits = pow(2, 10), M = 100;
    size_t N = 100000;

    Pulse rc_pulse(RC_function, {{"alpha", alpha}}, "RC");

    valarray<double> rc = rc_pulse.get_array(factor, fs);
    valarray<double> rc_fft = rc_pulse.get_array_fft(factor, fs);

    valarray<double> rc_fft_abs = abs(rc_fft)/fs;
    valarray<double> rc_fft_dB = 20*log10(rc_fft/fs);

    string filename = "rc_filter_t.txt";
    write_to_file(rc, filename);
    
    filename = "rc_filter_f.txt";
    write_to_file(rc_fft_dB, filename);

    valarray<double> eye_diagram_test = eye_diagram(N, rc, fs);

    filename = "rc_filter_eye_diagram.txt";
    write_to_file(eye_diagram_test, filename);

    valarray<double> rc_ber_isi = ber_isi(rc_pulse, fs, snr_db, nbits, M);

    printf("BER SSI ONLY\n");
    for (int i=0; i<rc_ber_isi.size(); i++)
    {
        printf("t/T = %.2f -> %.5e\n", OFFSETS[i], rc_ber_isi[i]);
    }

    valarray<double> rc_ber_cci = ber_cci(rc_pulse, fs, snr_db, sir_db, L, nbits, M);
    
    printf("BER CCI ONLY\n");
    for (int i=0; i<rc_ber_cci.size(); i++)
    {
        printf("t/T = %.2f -> %.5e\n", OFFSETS[i], rc_ber_cci[i]);
    }
}

// -- END OF FILE -- //