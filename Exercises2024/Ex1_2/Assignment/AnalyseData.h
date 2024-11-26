#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#pragma once


std::ifstream read_file(std::string input_file="./../input2D_float.txt");

std::vector<std::pair<double, double>> read_data_from_file(std::ifstream&);

void print_data(std::vector<std::pair<double, double>>);