#include <cmath>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <string>

#include "questions.hpp"
using namespace std;

int main() {  
    // double alpha = 0.35, omega = 0.10;
    // int snr_db = 15, sir_db = 10, L = 2, nbits = pow(2, 10), M = 100;
    // Define pulses and parameters
    const double mu = 1.6, epsilon = 0.1, gamma = 1, beta = 0.1;
    vector<function<double(double, unordered_map<string, double>)>> pulses = {RC_function, BTRC_function, IPLCP_function, ELP_function};
    vector<unordered_map<string, double>> params = {
        {},
        {},
        { {"mu", mu}, {"epsilon", epsilon}, {"gamma", gamma}},
        { {"beta", beta}}
    };
    vector<string> names = {"RC", "BTRC", "IPLCP", "ELP"};
    vector<double> alphas = {0.22, 0.5};
    vector<int> Ls = {2, 6};
    vector<int> snr_db = {10, 15};
    vector<int> sir_db = {10, 20};
    int factor = 4, fs = 10, NFFT = 1024;
    int L = 2, nbits = pow(2, 10), M = 100;
    int sir = 15;
    size_t N = 100000;
    // P1
    tuple<vector<valarray<double>>, vector<string>> P1_results = P1(pulses, params, names, alphas, factor, fs, NFFT);
    // Save results to CSV
    write_to_csv(get<0>(P1_results), get<1>(P1_results), "P1_results.csv");
    // P2
    tuple<vector<valarray<double>>, vector<string>, vector<double>> P2_results = P2(pulses, params, names, alphas[0], N, fs, factor);
    // Save results to CSV
    write_to_csv(get<0>(P2_results), get<1>(P2_results), "P2_results.csv");
    // Also save max distortion due to ISI
    // transform to valarray<double>
    valarray<double> max_distortion(get<2>(P2_results).data(), get<2>(P2_results).size());
    write_to_csv({max_distortion}, {"Max Distortion"}, "P3_results.csv");
    // P4
    // Define vectors for snr and alpha
    alphas = {0.22, 0.35, 0.5};
    tuple<vector<valarray<double>>, vector<string>> P4_results = P4(pulses, params, names, alphas, fs, snr_db, nbits, M);
    // Save results to CSV
    write_to_csv(get<0>(P4_results), get<1>(P4_results), "P4_results.csv");
    // P5
    tuple<vector<valarray<double>>, vector<string>> P5_results = P5(pulses, params, names, alphas, fs, snr_db[1], sir_db, Ls, nbits, M);
    // Save results to CSV
    write_to_csv(get<0>(P5_results), get<1>(P5_results), "P5_results.csv");
    // P6
    tuple<vector<valarray<double>>, vector<string>> P6_results = P6(pulses, params, names, alphas, fs, snr_db[1], sir, Ls[1], nbits, M);
    // Save results to CSV
    write_to_csv(get<0>(P6_results), get<1>(P6_results), "P6_results.csv");

    // P7
    vector<double> truncations = {5.0, 10.0};
    vector<int> sir_db2 = {10};
    vector<int> sir_db3 = {15};
    snr_db = {15};
    for (double truncation : truncations)
    {
        // Format truncation string
        stringstream truncation_ss;
        truncation_ss << fixed << setprecision(0) << truncation;
        // P7a
        tuple<vector<valarray<double>>, vector<string>> P7a_results = P4(pulses, params, names, alphas, fs, snr_db, nbits, M, truncation);
        // Save results to CSV
        write_to_csv(get<0>(P7a_results), get<1>(P7a_results), "P7a_results_truncation_" + truncation_ss.str() + ".csv");
        // P7b
        tuple<vector<valarray<double>>, vector<string>> P7b_results = P5(pulses, params, names, alphas, fs, snr_db[0], sir_db2, Ls, nbits, M, truncation);
        // Save results to CSV
        write_to_csv(get<0>(P7b_results), get<1>(P7b_results), "P7b_results_truncation_" + truncation_ss.str() + ".csv");
        // P7c
        tuple<vector<valarray<double>>, vector<string>> P7c_results = P6(pulses, params, names, {alphas[0]}, fs, snr_db[0], sir_db3[0], Ls[1], nbits, M, truncation);
        // Save results to CSV
        write_to_csv(get<0>(P7c_results), get<1>(P7c_results), "P7c_results_truncation_" + truncation_ss.str() + ".csv");
    }

}

// -- END OF FILE -- //