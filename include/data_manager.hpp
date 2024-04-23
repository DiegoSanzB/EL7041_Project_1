#ifndef DATA_MANAGER
#define DATA_MANAGER

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

#include "fft.hpp"

using namespace std;

const filesystem::path mainFilePath = filesystem::current_path();
const filesystem::path RESULTS_FOLDER_PATH = mainFilePath.parent_path() / "results";

void results_folder_check();
void write_to_file(CArray data, string filename);

#endif
// -- END OF FILE -- //