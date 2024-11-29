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

    bool go = true;
    while (go){
        std::cout<<"Please specify the function you want to run, enter 1 for print the data, enter 2 for calculate the magnitude , enter 3 for calculate x to the y, enter 4 for linear fitting, enter other int number to quit."<<std::endl;
        int command;
        std::cin >> command;
        switch (command)
        {
        case 1:{
            std::cout<<"How many lines would you like to print?";
            int n_lines;
            std::cin >> n_lines; 
            myprint(data,n_lines);
            break;
        }

        case 2:{
            std::vector<double> data_mag;
            //calculate the magnitude
            data_mag = mag(data);
            myprint(data_mag);
            my_save(data_mag, "mag.txt");
            std::cout<<"Results saved to mag.txt"<<std::endl;
            break;
        }

        case 3:{
            std::vector<double> x_power_y;
            x_power_y = x_to_the_y(data);
            myprint(x_power_y);
            my_save(x_power_y,"x_to_the_y.txt");
            std::cout<<"Results saved to x_to_the_y.txt"<<std::endl;
            break;
        }

        case 4:{
            std::pair<double, double> least_square_fit_result;
            least_square_fit_result = linearFit(data);
            double chi_square;
            chi_square = chi_squre_of_the_fit(least_square_fit_result,data);

            std::ofstream outFile("least_square_fit_result.txt");
            outFile << "y=" << least_square_fit_result.first << "x+" << least_square_fit_result.second;
            outFile.close();
            std::cout<<"Fit result is:" << "y=" << least_square_fit_result.first << "x+" << least_square_fit_result.second<<". Chi square of the fit is " << chi_square<<". Saved to least_square_fit_result.txt" <<std::endl;
            my_save(least_square_fit_result,chi_square,"fit_result.txt");
            std::cout<<"Results saved to fit_result.txt"<<std::endl;
            break;
        }
        
        default:{
            std::cout<<"Exiting."<<std::endl;
            go = false;
            break;
        }
        }// end switch

        
    }

    return 0;
}

