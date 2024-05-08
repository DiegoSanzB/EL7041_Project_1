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
    const double mu = 1.6, epsilon = 0.1, gamma = 1, beta = 0.1, omega = 0.10;
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
    int factor = 30, fs = 20, NFFT = pow(2, 14);
    int L = 2, nbits = pow(2, 10), M = 100;
    int sir = 15;
    size_t N = 100000;
    
    // P1
    tuple<vector<valarray<double>>, vector<valarray<double>>, vector<string>> P1_results = P1(pulses, params, names, alphas, factor, fs, NFFT);
    // Save results to CSV
    write_to_csv(get<0>(P1_results), get<2>(P1_results), "P1_time_results.csv");
    write_to_csv(get<1>(P1_results), get<2>(P1_results), "P1_freq_results.csv");
    printf("P1 computed!\n");
    
    // P2
    tuple<vector<valarray<double>>, vector<string>, vector<double>> P2_results = P2(pulses, params, names, alphas[0], N, fs, factor);
    // Save results to CSV
    write_to_csv(get<0>(P2_results), get<1>(P2_results), "P2_results.csv");
    printf("P2 computed!\n");
    
    // P3
    // Also save max distortion due to ISI
    // transform to valarray<double>
    valarray<double> max_distortion(get<2>(P2_results).data(), get<2>(P2_results).size());
    write_to_csv({max_distortion}, {"Max Distortion"}, "P3_results.csv");
    printf("P3 computed!\n");
    
    // P4
    // Define vectors for snr and alpha
    alphas = {0.22, 0.35, 0.5};
    tuple<vector<valarray<double>>, vector<string>> P4_results = P4(pulses, params, names, alphas, omega, snr_db, nbits, M);
    // Save results to CSV
    write_to_csv(get<0>(P4_results), get<1>(P4_results), "P4_results.csv");
    printf("P4 computed!\n");

    // P5
    tuple<vector<valarray<double>>, vector<string>> P5_results = P5(pulses, params, names, alphas, omega, snr_db[1], sir_db, Ls, nbits, M);
    // Save results to CSV
    write_to_csv(get<0>(P5_results), get<1>(P5_results), "P5_results.csv");
    printf("P5 computed!\n");
    
    // P6
    tuple<vector<valarray<double>>, vector<string>> P6_results = P6(pulses, params, names, alphas, omega, snr_db[1], sir, Ls[1], nbits, M);
    // Save results to CSV
    write_to_csv(get<0>(P6_results), get<1>(P6_results), "P6_results.csv");
    printf("P6 computed!\n");

    // P7
    vector<double> truncations = {5.0, 10.0};
    vector<int> sir_db2 = {10};
    vector<int> sir_db3 = {15};
    snr_db = {10};
    for (double truncation : truncations)
    {
        // Format truncation string
        stringstream truncation_ss;
        truncation_ss << fixed << setprecision(0) << truncation;
        // P7a
        tuple<vector<valarray<double>>, vector<string>> P7a_results = P4(pulses, params, names, alphas, omega, snr_db, nbits, M, truncation);
        // Save results to CSV
        write_to_csv(get<0>(P7a_results), get<1>(P7a_results), "P7a_results_truncation_" + truncation_ss.str() + ".csv");
        // P7b
        tuple<vector<valarray<double>>, vector<string>> P7b_results = P5(pulses, params, names, alphas, omega, sir_db3[0], sir_db2, Ls, nbits, M, truncation);
        // Save results to CSV
        write_to_csv(get<0>(P7b_results), get<1>(P7b_results), "P7b_results_truncation_" + truncation_ss.str() + ".csv");
        // P7c
        tuple<vector<valarray<double>>, vector<string>> P7c_results = P6(pulses, params, names, {alphas[0]}, omega, sir_db3[0], sir_db3[0], Ls[1], nbits, M, truncation);
        // Save results to CSV
        write_to_csv(get<0>(P7c_results), get<1>(P7c_results), "P7c_results_truncation_" + truncation_ss.str() + ".csv");
    }
    printf("P7 computed!\n");

    return 0;
}

// -- END OF FILE -- //