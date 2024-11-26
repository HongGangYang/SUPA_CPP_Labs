#include "AnalyseData.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>


std::ifstream read_file(std::string input_file){
    std::ifstream data_file(input_file);
    
    // Check that the file is opened properly
    if (!data_file.is_open()) {
        std::cout << "Failed to open the file." << std::endl;
        // return 1; 
    }

    return data_file;
}

std::vector<std::pair<double, double>> read_data_from_file(std::ifstream& data_file){
    std::vector<std::pair<double, double>> data;
    std::string line;

    //to skip the first line which is "x,y" rather than the data.
    std::getline(data_file, line);

    while (std::getline(data_file, line)) {
        std::stringstream ss(line);
        std::string value;
        double x, y;

        // Read x
        std::getline(ss, value, ',');
        x = std::stod(value);

        // Read y
        std::getline(ss, value, ',');
        y = std::stod(value);

        data.emplace_back(x, y); // Add the pair to the vector
    }
    return data;
}

void print_data(std::vector<std::pair<double, double>> data){
    for (const auto& point : data) {
        std::cout <<  point.first << "," << point.second << std::endl;
    }
}