// Created by Hong-Gang on 11th Dec for the assignment
#include "FiniteFunctions.h"
#include <iostream>
#include <cmath>
#include <stdexcept>

int main(){ 
    // Gaussian test_;
    // test_ = Gaussian(-6.,6.,2,0.85,"testGaussian");
    // Cauchy_Lorentz test_;
    // test_ = Cauchy_Lorentz(-6.,6.,2,0.5,"testCauchyLorentz");
    Negative_Crystal_Ball test_;
    test_ = Negative_Crystal_Ball(-6.,6.,2,3,0.469,1,"testNegativeCrystalBall");
    std::vector<double> random_data;
    random_data = readDataFromFile();
    test_.plotData(random_data,100,true);
    test_.plotFunction();
    return 0;
}