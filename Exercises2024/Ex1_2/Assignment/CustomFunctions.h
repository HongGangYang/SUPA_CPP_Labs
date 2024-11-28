#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#pragma once


std::ifstream read_file(std::string input_file="./../input2D_float.txt");

std::vector<std::pair<double, double>> read_data_from_file(std::ifstream&);

std::vector<double> mag(std::vector<std::pair<double, double>>);

void myprint(std::vector<std::pair<double, double>>,int);

void myprint(std::vector<double>);

std::pair<double, double> linearFit(const std::vector<std::pair<double, double>>&);

double chi_squre_of_the_fit(std::pair<double, double>, std::vector<std::pair<double, double>>, std::string file_path_for_sigma = "./../error2D_float.txt");

double x_to_the_y_single_point(double, int);

std::vector<double> x_to_the_y(std::vector<std::pair<double, double>>);