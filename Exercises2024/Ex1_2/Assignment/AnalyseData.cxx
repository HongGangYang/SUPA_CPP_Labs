// Created by Hong-Gang on 26th Nov for the assignment
#include "AnalyseData.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

int main() {

    std::ifstream data_file = read_file();

    std::vector<std::pair<double, double>> data;

    data = read_data_from_file(data_file);

    print_data(data);

    return 0;
}

