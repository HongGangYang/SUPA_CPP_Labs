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
    std::cout<<"Reading data finished."<<std::endl;

    //test
    std::vector<double> x_power_y;
    x_power_y = x_to_the_y(data);
    myprint(x_power_y);


    // declare the magnitude vector
    std::vector<double> data_mag;
    //calculate the magnitude
    data_mag = mag(data);
    std::cout<<"Magnitude calculation finished."<<std::endl;

    std::pair<double, double> least_square_fit_result;
    least_square_fit_result = linearFit(data);
    double chi_square;
    chi_square = chi_squre_of_the_fit(least_square_fit_result,data);

    std::ofstream outFile("least_square_fit_result.txt");
    outFile << "y=" << least_square_fit_result.first << "x+" << least_square_fit_result.second;
    outFile.close();
    std::cout<<"Fit result is:" << "y=" << least_square_fit_result.first << "x+" << least_square_fit_result.second<<". Chi square of the fit is " << chi_square<<". Saved to least_square_fit_result.txt" <<std::endl;

    while (true){
        std::cout<<"Please specify the function you want to run, choose among \"print_data\",\"print_mag\", and \"quit\"."<<std::endl;
        std::string command;
        std::cin >> command;
        if (command == "print_data"){
            std::cout<<"How many lines would you like to print?";
            int n_lines;
            std::cin >> n_lines; 
            myprint(data,n_lines);
        }
        else if (command == "print_mag"){
            myprint(data_mag);
        }
        else if (command == "quit"){
            break;
        }
        else {
            std::cout<<"command not valid"<<std::endl;
        }
    }


    // // print the data
    // print_data(data,100);
    // // print the magniture
    // print_mag(data_mag);

    return 0;
}

