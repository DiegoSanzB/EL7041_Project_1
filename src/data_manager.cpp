#include "data_manager.hpp"

void results_folder_check()
{
    if (!filesystem::exists(RESULTS_FOLDER_PATH))
    {
        filesystem::create_directories(RESULTS_FOLDER_PATH);
    }
}

void write_to_file(valarray<double> data, string filename)
{
    results_folder_check();
    ofstream output_file;
    output_file.open(RESULTS_FOLDER_PATH / filename);

    for (Complex element : data)
    {
        output_file << element << endl;
    }
}