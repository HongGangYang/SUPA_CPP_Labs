#include "CustomFunctions.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
 

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

std::vector<double> mag(std::vector<std::pair<double, double>> data){
    std::vector<double> mag_;
    for (const auto& point : data) {
        double mag_temp;
        mag_temp = pow(point.first*point.first+point.second*point.second,0.5);
        mag_.push_back(mag_temp);
    }
    return mag_;
}

void myprint(std::vector<std::pair<double, double>> data, int n){
    if (n > data.size()){
        std::cout<<"Warning: n is bigger than the size of the data. The first 5 data is printed instead."<<std::endl;
        std::vector<std::pair<double, double>> first_n(data.begin(), data.begin() + 5);
        for (const auto& point : first_n) {
         std::cout <<  point.first << "," << point.second << std::endl;
    } 
    }
    else {
        std::vector<std::pair<double, double>> first_n(data.begin(), data.begin() + n);
        for (const auto& point : first_n) {
        std::cout <<  point.first << "," << point.second << std::endl;
    }
    }
    
}

void myprint(std::vector<double> mag){
    for (auto& mag_ : mag) {
    std::cout <<  mag_ << std::endl;
    }
    
}

std::pair<double, double> linearFit(const std::vector<std::pair<double, double>>& data) {
    double sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_xy = 0.0;
    int N = data.size();
    
    // Compute the sums
    for (const auto& point : data) {
        double x = point.first;
        double y = point.second;

        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
    }

    // Calculate slope (p) and intercept (q)
    double denominator = N * sum_xx - sum_x * sum_x;

    double p = (N * sum_xy - sum_x * sum_y) / denominator;
    double q = (sum_xx*sum_y-sum_xy*sum_x) / denominator;

    return {p, q}; 
}

double chi_squre_of_the_fit(std::pair<double, double> fit_result, std::vector<std::pair<double, double>> data, std::string file_path_for_sigma){
    std::ifstream sigma_file = read_file(file_path_for_sigma);
    std::vector<std::pair<double, double>> sigma;
    sigma = read_data_from_file(sigma_file);

    int n_points = data.size();
    double chi_square = 0.;
    for (int i = 0; i<n_points ; i++) {
        chi_square += pow((data[i].second-fit_result.first*data[i].first-fit_result.second),2)/pow(sigma[i].second,2);
        }
    chi_square /= (n_points-2);
    return chi_square;
}

double x_to_the_y_single_point(double x, int y){
    if (y==0){
        return 1.0;
    }
    return x*x_to_the_y_single_point(x,y-1);
}

std::vector<double> x_to_the_y(std::vector<std::pair<double, double>> data){
    int n_points = data.size();
    std::vector<double> result;
    for (int i = 0; i<n_points ; i++) {
        result.push_back(x_to_the_y_single_point(data[i].first,std::round(data[i].second)));
    }
    return result;
}



void my_save(const std::vector<double>& data, const std::string& filename) {
    std::ofstream output_file(filename); 
    for (const double& value : data) {
        output_file << value << "\n"; // Write each element on a new line
    }

    output_file.close(); // Close the file
}

void my_save(std::pair<double, double>& fit_result_, double chi_square_, const std::string& filename){
    std::ofstream output_file(filename); 
    output_file<<"Fit result is:" << "y=" << fit_result_.first << "x+" << fit_result_.second<<". Chi square of the fit is " << chi_square_<<std::endl;
    output_file.close(); // Close the file
}