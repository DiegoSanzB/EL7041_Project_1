#ifndef DATA_MANAGER
#define DATA_MANAGER

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <valarray>


// #include "fft.hpp"

using namespace std;
const filesystem::path mainFilePath = filesystem::current_path();
const filesystem::path RESULTS_FOLDER_PATH = mainFilePath.parent_path() / "results"; // Update with your folder path

void write_to_csv(const vector<valarray<double>> &data, const vector<string> &column_names, string CSV_FILENAME);

#endif
// -- END OF FILE -- //