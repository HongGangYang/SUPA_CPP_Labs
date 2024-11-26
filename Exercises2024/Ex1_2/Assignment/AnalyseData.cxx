// Created by Hong-Gang on 26th Nov for the assignment
// g++ -std=c++20 -w CustomFunctions.cxx AnalyseData.cxx -o test
#include "CustomFunctions.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

int main() {
    // open the file
    std::ifstream data_file = read_file();
    // declare the data read from the file
    std::vector<std::pair<double, double>> data;
    // read the data from the file
    data = read_data_from_file(data_file);
    // declare the magnitude vector
    std::vector<double> data_mag;
    //calculate the magnitude
    data_mag = mag(data);
    // print the data
    print_data(data,100);
    // print the magniture
    print_mag(data_mag);

    return 0;
}

