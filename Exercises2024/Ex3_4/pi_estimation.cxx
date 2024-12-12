// Created by Hong-Gang on 12th Dec for the assignment
// g++ -std=c++20 -o pi_estimation pi_estimation.cxx
#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>

double estimate_pi(double radius, int n_points);

int main(){
    std::cout<<"Please enter the radius:";
    double radius_;
    std::cin >> radius_;
    std::cout<<"Please enter the random points you want to generate:";
    int n_points_;
    std::cin >> n_points_;
    std::cout<<"The estimation for pi given the setting is:"<<std::fixed<<std::setprecision(10)<<estimate_pi(radius_,n_points_)<<std::endl;
    return 0;
}

double estimate_pi(double radius, int n_points) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> uniform_dist(-radius, radius);

    int n_points_indide = 0;

    for (int i = 0; i < n_points; ++i) {
        double x = uniform_dist(gen);
        double y = uniform_dist(gen);

        if (x * x + y * y <= radius * radius) {
            n_points_indide++;
        }
    }

    double pi = 4.0 * n_points_indide / n_points;
    return pi;
}