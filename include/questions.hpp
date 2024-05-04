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
#include "impulses.hpp"

// #ifndef QUESTIONS
// #define QUESTIONS

using namespace std;

tuple<vector<valarray<double>>, vector<string>> P1(vector<function<double(double, unordered_map<string, double>)>> pulses, vector<unordered_map<string, double>> params, vector<string> names, vector<double> alphas, int factor, int fs, int NFFT);
tuple<vector<valarray<double>>, vector<string>, vector<double>> P2(vector<function<double(double, unordered_map<string, double>)>> pulses, vector<unordered_map<string, double>> params, vector<string> names, double alpha, size_t N, int fs, int factor);
tuple<vector<valarray<double>>, vector<string>> P4(vector<function<double(double, unordered_map<string, double>)>> pulses, vector<unordered_map<string, double>> params, vector<string> names, vector<double> alphas, int fs, vector<int> snr_db, int nbits, int M, double truncation = numeric_limits<double>::infinity());
tuple<vector<valarray<double>>, vector<string>> P5(vector<function<double(double, unordered_map<string, double>)>> pulses, vector<unordered_map<string, double>> params, vector<string> names, vector<double> alphas, int fs, int snr_db, vector<int> sir_db, vector<int> L, int nbits, int M, double truncation = numeric_limits<double>::infinity());
// tuple<vector<valarray<double>>, vector<string>> P6(vector<function<double(double, unordered_map<string, double>)>> pulses, vector<unordered_map<string, double>> params, vector<string> names, vector<double> alphas, int fs, int snr_db, int sir_db, int L, int nbits, int M, double truncation);
// #endif

// -- END OF FILE -- //