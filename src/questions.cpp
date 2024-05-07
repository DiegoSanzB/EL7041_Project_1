#include "questions.hpp"

// P1 - Plot pulses in time and frequency domain for different alpha values and pulse shapes
// Inputs:
// - pulses: vector of pulse functions
// - params: vector of pulse parameters for each pulse
// - names: vector of pulse names
// - alphas: vector of alpha values
// - factor: factor to determine the array size
// - fs: sampling frequency
// Outputs:
// - vector of valarray<double> containing the pulses in time domain and frequency domain
// - vector of strings containing the pulse names
// Note: in I give n pulses, I should have n*len(alphas)*2*(valarray<double>) in the output
tuple<vector<valarray<double>>, vector<valarray<double>>, vector<string>> P1(vector<function<double(double, unordered_map<string, double>)>> pulses, vector<unordered_map<string, double>> params, vector<string> names, vector<double> alphas, int factor, int fs, int NFFT){
    // Create an empty vector to store the results
    vector<valarray<double>> time_results;
    vector<valarray<double>> freq_results;
    vector<string> pulse_names;
    for (int i = 0; i < pulses.size(); i++){
        for (int j = 0; j < alphas.size(); j++){
            // Include alpha in params
            params[i]["alpha"] = alphas[j];
            // create a string alpha with specific precision
            stringstream ss;
            ss << fixed << setprecision(2) << alphas[j];
            // Create a name for the pulse
            string pulse_name = names[i] + "_alpha_" + ss.str();
            pulse_names.push_back(pulse_name);
            // Create pulse object
            Pulse pulse(pulses[i], params[i], pulse_name);
            // Get pulse array in time domain
            valarray<double> pulse_time = pulse.get_array(factor, fs);
            // Get pulse array in frequency domain
            valarray<double> pulse_freq = pulse.get_array_fft(factor, fs, NFFT);
            // Store the results
            time_results.push_back(pulse_time);
            freq_results.push_back(pulse_freq);
        } 
    }
    return make_tuple(time_results, freq_results, pulse_names);
}

// P2 - Plot eye diagram for the pulses and a fixed alpha value
// Inputs:
// - pulses: vector of pulse functions
// - params: vector of pulse parameters for each pulse
// - names: vector of pulse names
// - alpha: alpha value
// - N: number of symbols
// - fs: sampling frequency
// - factor: factor to determine the array size
// Outputs:
// - vector of valarray<double> containing the eye diagrams
// - vector of strings containing the pulse names
// Also, return the max distortion due to ISI so P3 is done
tuple<vector<valarray<double>>, vector<string>, vector<double>> P2(vector<function<double(double, unordered_map<string, double>)>> pulses, vector<unordered_map<string, double>> params, vector<string> names, double alpha, size_t N, int fs, int factor){
    // Create an empty vector to store the results
    vector<valarray<double>> results;
    vector<string> pulse_names;
    vector<double> max_distortions;
    // create a string alpha with specific precision
    stringstream ss;
    ss << fixed << setprecision(2) << alpha;
    for (int i = 0; i < pulses.size(); i++){
        // Include alpha in params
        params[i]["alpha"] = alpha;
        // Create a name for the pulse
        string pulse_name = names[i] + "_alpha_" + ss.str();
        pulse_names.push_back(pulse_name);
        // Create pulse object
        Pulse pulse(pulses[i], params[i], pulse_name);
        // Create pulse array in time domain
        valarray<double> pulse_array = pulse.get_array(factor, fs);
        // Get eye diagram
        valarray<double> eye_diagram_values = eye_diagram(N, pulse_array, fs);
        // Store the results
        results.push_back(eye_diagram_values);
        // Get max distortion due to ISI
        double max_dist = abs(eye_diagram_values).max();
        max_distortions.push_back(max_dist);    
    }
    return make_tuple(results, pulse_names, max_distortions);
}

// Everything BER related 

// P4 - Plot BER for ISI only for the pulses and multiple alpha values
// Inputs:
// - pulses: vector of pulse functions
// - params: vector of pulse parameters for each pulse
// - names: vector of pulse names
// - alphas: vector of alpha values
// - fs: sampling frequency
// - snr_db:  vector of SNR values in dB
// - nbits: number of bits
// - M: not sure ask uwu
// truncation: truncation value for the pulse
// Outputs:
// - vector of valarray<double> containing the BER values for offset values
// - vector of strings containing the pulse names
tuple<vector<valarray<double>>, vector<string>> P4(vector<function<double(double, unordered_map<string, double>)>> pulses, vector<unordered_map<string, double>> params, vector<string> names, vector<double> alphas, double omega, vector<int> snr_db, int nbits, int M, double truncation){
    // Create an empty vector to store the results
    vector<valarray<double>> results;
    vector<string> pulse_names;
    for (int i = 0; i < pulses.size(); i++){
        for (int j = 0; j < alphas.size(); j++){
            // Include alpha in params
            params[i]["alpha"] = alphas[j];
            // Create a string for alpha
            stringstream ss_alpha;
            ss_alpha << fixed << setprecision(2) << alphas[j];
            for (int k = 0; k < snr_db.size(); k++){
                // Create string for snr
                stringstream ss_snr;
                ss_snr << snr_db[k];    
                // Create a name for the pulse
                string pulse_name = names[i] + "_alpha_" + ss_alpha.str() + "_snr_" + ss_snr.str();
                pulse_names.push_back(pulse_name);
                // Create pulse object
                Pulse pulse(pulses[i], params[i], pulse_name);
                // Get BER values for ISI
                valarray<double> ber_isi_result = ber_isi(pulse, omega, snr_db[k], nbits, M, truncation);
                // Store the results
                results.push_back(ber_isi_result);
            }
        } 
    }
    return make_tuple(results, pulse_names);
}
// P5 - Plot BER for CCI only for the pulses with multiple alpha values, L values and SIR values
// Inputs:
// - pulses: vector of pulse functions
// - params: vector of pulse parameters for each pulse
// - names: vector of pulse names
// - alphas: vector of alpha values
// - fs: sampling frequency
// - snr_db:  snr value in dB
// - sir_db: vector of SIR values in dB
// - L: vector of L values
// - nbits: number of bits
// - M: not sure ask uwu
// truncation: truncation value for the pulse
// Outputs:
// - vector of valarray<double> containing the BER values for offset values
// - vector of strings containing the pulse names
tuple<vector<valarray<double>>, vector<string>> P5(vector<function<double(double, unordered_map<string, double>)>> pulses, vector<unordered_map<string, double>> params, vector<string> names, vector<double> alphas, double omega, int snr_db, vector<int> sir_db, vector<int> L, int nbits, int M, double truncation){
    // Create an empty vector to store the results
    vector<valarray<double>> results;
    vector<string> pulse_names;
    for (int i = 0; i < pulses.size(); i++){
        for (int j = 0; j < alphas.size(); j++){
            // Include alpha in params
            params[i]["alpha"] = alphas[j];
            // Create a string for alpha
            stringstream ss_alpha;
            ss_alpha << fixed << setprecision(2) << alphas[j];
            for (int k = 0; k < sir_db.size(); k++){
                // Create a string for sir
                stringstream ss_sir;
                ss_sir << sir_db[k];
                for (int l = 0; l < L.size(); l++){
                    // Create a name for the pulse
                    string pulse_name = names[i] + "_alpha_" + ss_alpha.str() + "_sir_" + ss_sir.str() + "_L_" + to_string(L[l]);
                    pulse_names.push_back(pulse_name);
                    // Create pulse object
                    Pulse pulse(pulses[i], params[i], pulse_name);
                    // Get BER values for CCI
                    valarray<double> ber_cci_result = ber_cci(pulse, omega, snr_db, sir_db[k], L[l], nbits, M, truncation);
                    // Store the results
                    results.push_back(ber_cci_result);
                }
            }
        } 
    }
    return make_tuple(results, pulse_names);
}

// P6 - Plot BER for ISI and CCI for the pulses with multiple alpha values, fixed L value, SNR and SIR values
// Inputs:
// - pulses: vector of pulse functions
// - params: vector of pulse parameters for each pulse
// - names: vector of pulse names
// - alphas: vector of alpha values
// - fs: sampling frequency
// - snr_db:  SNR value in dB
// - sir_db: SIR value in dB
// - L: fixed L value
// - nbits: number of bits
// - M: not sure ask uwu
// Outputs:
// - vector of valarray<double> containing the BER values for offset values
// - vector of strings containing the pulse names
tuple<vector<valarray<double>>, vector<string>> P6(vector<function<double(double, unordered_map<string, double>)>> pulses, vector<unordered_map<string, double>> params, vector<string> names, vector<double> alphas,  double omega, int snr_db, int sir_db, int L, int nbits, int M, double truncation){
    // Create an empty vector to store the results
    vector<valarray<double>> results;
    vector<string> pulse_names;
    for (int i = 0; i < pulses.size(); i++){
        for (int j = 0; j < alphas.size(); j++){
            // Create a string for alpha
            stringstream ss_alpha;
            ss_alpha << fixed << setprecision(2) << alphas[j];
            // Include alpha in params
            params[i]["alpha"] = alphas[j];
            // Create a name for the pulse
            string pulse_name = names[i] + "_alpha_" + ss_alpha.str();
            pulse_names.push_back(pulse_name);
            // Create pulse object
            Pulse pulse(pulses[i], params[i], pulse_name);
            // Get BER values for ISI and CCI
            valarray<double> ber_isi_cci_result = ber_isi_cci(pulse, omega, snr_db, sir_db, L, nbits, M, truncation);
            // Store the results
            results.push_back(ber_isi_cci_result);
        }
    }
    return make_tuple(results, pulse_names);
}



// P7 - Same process with truncated pulses, many truncation limits
// Order it in main c:


