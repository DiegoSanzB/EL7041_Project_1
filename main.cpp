#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>
#include <string>

#include "data_manager.hpp"
#include "eye_diagram.hpp"
#include "fft.hpp"
#include "impulses/impulses.hpp"

using namespace std;

int main() {  
    float alpha = 0.35;

    float fs = 10;
    int NFFT = 1024;
    size_t N = 100000;

    CArray rc_filter = raised_cosine(alpha, fs, NFFT);

    CArray rc_filter_spectrum = fft(rc_filter);

    CArray rc_fft = abs(rc_filter_spectrum)/fs;
    CArray rc_fft_dB = 20*log10(rc_filter_spectrum/fs);

    string filename = "rc_filter_t.txt";
    write_to_file(rc_filter, filename);
    
    filename = "rc_filter_f.txt";
    write_to_file(rc_fft, filename);

    CArray eye_diagram_test = eye_diagram(N, rc_filter, fs);

    filename = "rc_filter_eye_diagram.txt";
    write_to_file(eye_diagram_test, filename);
}

