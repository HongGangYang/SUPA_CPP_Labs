#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include "gnuplot-iostream.h"

#pragma once //Replacement for IFNDEF

class FiniteFunction{   

public:
  FiniteFunction(); //Empty constructor
  FiniteFunction(double range_min, double range_max, std::string outfile); //Variable constructor
  ~FiniteFunction(); //Destructor
  double rangeMin(); //Low end of the range the function is defined within
  double rangeMax(); //High end of the range the function is defined within
  double integral(int Ndiv = 1000); 
  std::vector< std::pair<double,double> > scanFunction(int Nscan = 1000); //Scan over function to plot it (slight hack needed to plot function in gnuplot)
  void setRangeMin(double RMin);
  void setRangeMax(double RMax);
  void setOutfile(std::string outfile);
  void plotFunction(); //Plot the function using scanFunction
  
  //Plot the supplied data points (either provided data or points sampled from function) as a histogram using NBins
  void plotData(std::vector<double> &points, int NBins, bool isdata=true); //NB! use isdata flag to pick between data and sampled distributions
  virtual void printInfo(); //Dump parameter info about the current function (Overridable)
  virtual double callFunction(double x); //Call the function with value x (Overridable)

  //Protected members can be accessed by child classes but not users
protected:
  double m_RMin;
  double m_RMax;
  double m_Integral;
  int m_IntDiv = 0; //Number of division for performing integral
  std::string m_FunctionName;
  std::string m_OutData; //Output filename for data
  std::string m_OutPng; //Output filename for plot
  std::vector< std::pair<double,double> > m_data; //input data points to plot
  std::vector< std::pair<double,double> > m_samples; //Holder for randomly sampled data 
  std::vector< std::pair<double,double> > m_function_scan; //holder for data from scanFunction (slight hack needed to plot function in gnuplot)
  bool m_plotfunction = false; //Flag to determine whether to plot function
  bool m_plotdatapoints = false; //Flag to determine whether to plot input data
  bool m_plotsamplepoints = false; //Flag to determine whether to plot sampled data 
  double integrate(int Ndiv);
  std::vector< std::pair<double, double> > makeHist(std::vector<double> &points, int Nbins); //Helper function to turn data points into histogram with Nbins
  void checkPath(std::string outstring); //Helper function to ensure data and png paths are correct
  void generatePlot(Gnuplot &gp); 
  
private:
  double invxsquared(double x); //The default functional form
};

std::vector<double> readDataFromFile(const std::string& filePath="./Outputs/data/MysteryData15041.txt");

class Gaussian : public FiniteFunction {
public:
  Gaussian();
  Gaussian(double range_min, double range_max,double mu_ ,double sigma_,std::string outfile);
  virtual double callFunction(double x); //Call the function with value x (Overridable)

protected:
  double m_mu;
  double m_sigma;

private:
  double gaussian_function(double x, double mu, double sigma);
};

class Cauchy_Lorentz : public FiniteFunction {
public:
  Cauchy_Lorentz();
  Cauchy_Lorentz(double range_min, double range_max,double x_0 ,double gamma_,std::string outfile);
  virtual double callFunction(double x); //Call the function with value x (Overridable)

protected:
  double m_x_0;
  double m_gamma;

private:
  double Cauchy_Lorentz_function(double x, double x_0, double gamma);
};

class Negative_Crystal_Ball : public FiniteFunction {
public:
  Negative_Crystal_Ball();
  Negative_Crystal_Ball(double range_min, double range_max, double x_bar_, double n_, double sigma_, double alpha_, std::string outfile);
  virtual double callFunction(double x); //Call the function with value x (Overridable)

protected:
  double m_x_bar;
  double m_n;
  double m_sigma;
  double m_alpha;
  double m_A;
  double m_B;
  double m_C;
  double m_D;
  double m_N;

private:
  double Negative_Crystal_Ball_function(double x, double x_bar_, double n_, double sigma_, double alpha_, double A_, double B_, double N_);
};