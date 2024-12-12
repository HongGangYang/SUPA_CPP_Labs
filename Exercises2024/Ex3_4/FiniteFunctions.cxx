#include <iostream>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include "FiniteFunctions.h"
#include <filesystem> //To check extensions in a nice way

#include "gnuplot-iostream.h" //Needed to produce plots (not part of the course) 

using std::filesystem::path;

//Empty constructor
FiniteFunction::FiniteFunction(){
  m_RMin = -5.0; 
  m_RMax = 5.0;
  this->checkPath("DefaultFunction");
  m_Integral = NULL;
}

Gaussian::Gaussian(){
  m_RMin = -5.0; 
  m_RMax = 5.0;
  m_mu = 0;
  m_sigma = 1;
  this->checkPath("DefaultFunction");
  m_Integral = NULL;
}

Cauchy_Lorentz::Cauchy_Lorentz(){
  m_RMin = -5.0; 
  m_RMax = 5.0;
  m_x_0 = 2;
  m_gamma = 1;
  this->checkPath("DefaultFunction");
  m_Integral = NULL;
}

Negative_Crystal_Ball::Negative_Crystal_Ball(){
  m_RMin = -5.0; 
  m_RMax = 5.0;
  m_x_bar = 2;
  m_n = 2;
  m_sigma = 1;
  m_alpha = 2;
  m_A = std::pow(m_n / std::fabs(m_alpha), m_n) * std::exp(-0.5 * m_alpha * m_alpha);
  m_B = m_n / std::fabs(m_alpha) - std::fabs(m_alpha);
  m_C = m_n / std::fabs(m_alpha) * (1.0 / (m_n - 1.0)) * std::exp(-0.5 * m_alpha * m_alpha);
  m_D = std::sqrt(M_PI / 2.0) * (1 + std::erf(std::fabs(m_alpha) / std::sqrt(2.0)));
  m_N = 1.0 / (m_sigma * (m_C + m_D));
  this->checkPath("DefaultFunction");
  m_Integral = NULL;
}

//initialised constructor
FiniteFunction::FiniteFunction(double range_min, double range_max, std::string outfile){
  m_RMin = range_min;
  m_RMax = range_max;
  m_Integral = NULL;
  this->checkPath(outfile); //Use provided string to name output files
}

Gaussian::Gaussian(double range_min, double range_max,double mu_ ,double sigma_, std::string outfile){
  m_RMin = range_min;
  m_RMax = range_max;
  m_mu = mu_;
  m_sigma = sigma_;
  m_Integral = NULL;
  this->checkPath(outfile); //Use provided string to name output files
}

Cauchy_Lorentz::Cauchy_Lorentz(double range_min, double range_max,double x_0_ ,double gamma_, std::string outfile){
  m_RMin = range_min;
  m_RMax = range_max;
  m_x_0 = x_0_;
  m_gamma = gamma_;
  m_Integral = NULL;
  this->checkPath(outfile); //Use provided string to name output files
}

Negative_Crystal_Ball::Negative_Crystal_Ball(double range_min, double range_max, double x_bar_, double n_, double sigma_, double alpha_, std::string outfile){
  m_RMin = range_min;
  m_RMax = range_max;
  m_x_bar = x_bar_;
  m_n = n_;
  m_sigma = sigma_;
  m_alpha = alpha_;
  m_A = std::pow(m_n / std::fabs(m_alpha), m_n) * std::exp(-0.5 * m_alpha * m_alpha);
  m_B = m_n / std::fabs(m_alpha) - std::fabs(m_alpha);
  m_C = m_n / std::fabs(m_alpha) * (1.0 / (m_n - 1.0)) * std::exp(-0.5 * m_alpha * m_alpha);
  m_D = std::sqrt(M_PI / 2.0) * (1 + std::erf(std::fabs(m_alpha) / std::sqrt(2.0)));
  m_N = 1.0 / (m_sigma * (m_C + m_D));
  m_Integral = NULL;
  this->checkPath(outfile); //Use provided string to name output files
}

//Plots are called in the destructor
//SUPACPP note: They syntax of the plotting code is not part of the course
FiniteFunction::~FiniteFunction(){
  Gnuplot gp; //Set up gnuplot object
  this->generatePlot(gp); //Generate the plot and save it to a png using "outfile" for naming 
}

/*
###################
//Setters
###################
*/ 
void FiniteFunction::setRangeMin(double RMin) {m_RMin = RMin;};
void FiniteFunction::setRangeMax(double RMax) {m_RMax = RMax;};
void FiniteFunction::setOutfile(std::string Outfile) {this->checkPath(Outfile);};

/*
###################
//Getters
###################
*/ 
double FiniteFunction::rangeMin() {return m_RMin;};
double FiniteFunction::rangeMax() {return m_RMax;};

/*
###################
//Function eval
###################
*/ 
double FiniteFunction::invxsquared(double x) {return 1/(1+x*x);};
double FiniteFunction::callFunction(double x) {return this->invxsquared(x);}; //(overridable)

double Gaussian::gaussian_function(double x, double mu, double sigma) {return pow(2.7183,-pow((x-mu)/sigma,2))/(sigma*pow(6.28,0.5));};
double Gaussian::callFunction(double x) {return this->gaussian_function(x,m_mu,m_sigma);};

double Cauchy_Lorentz::Cauchy_Lorentz_function(double x, double x_0, double gamma) {return 1/(3.14*gamma*(1+pow((x-x_0)/gamma,2)));};
double Cauchy_Lorentz::callFunction(double x) {return this->Cauchy_Lorentz_function(x,m_x_0,m_gamma);};

double Negative_Crystal_Ball::Negative_Crystal_Ball_function(double x, double x_bar_, double n_, double sigma_, double alpha_, double A_, double B_, double N_) {
  double t = (x - x_bar_) / sigma_;
  if (t > -alpha_) {
            return N_ * std::exp(-0.5 * t * t);
        } else {
            return N_ * A_ * std::pow(B_ - t, -n_);
        }
};
double Negative_Crystal_Ball::callFunction(double x) {return this->Negative_Crystal_Ball_function(x,m_x_bar,m_n,m_sigma,m_alpha,m_A,m_B,m_N);};

/*
###################
Integration by hand (output needed to normalise function when plotting)
###################
*/ 
double FiniteFunction::integrate(int Ndiv){ //private
  //ToDo write an integrator
  //function should be very close to zero at m_Rmin and m_Rmax
  double temp = 0.;
  double step = (m_RMax - m_RMin)/(double)Ndiv;
  for (int i = 0; i < Ndiv; i++) {
        double x1 = m_RMin + i * step;
        double x2 = x1 + step;
        double mid = (x1 + x2) / 2;
        temp += this->callFunction(mid) * step;
    }

  return temp;  
}
double FiniteFunction::integral(int Ndiv) { //public
  if (Ndiv <= 0){
    std::cout << "Invalid number of divisions for integral, setting Ndiv to 1000" <<std::endl;
    Ndiv = 1000;
  }
  if (m_Integral == NULL || Ndiv != m_IntDiv){
    m_IntDiv = Ndiv;
    m_Integral = this->integrate(Ndiv);
    return m_Integral;
  }
  else return m_Integral; //Don't bother re-calculating integral if Ndiv is the same as the last call
}

/*
###################
//Helper functions 
###################
*/
// Generate paths from user defined stem
void FiniteFunction::checkPath(std::string outfile){
 path fp = outfile;
 m_FunctionName = fp.stem(); 
 m_OutData = m_FunctionName+".data";
 m_OutPng = m_FunctionName+".png";
}

//Print (overridable)
void FiniteFunction::printInfo(){
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
}

/*
###################
//Plotting
###################
*/

//Hack because gnuplot-io can't read in custom functions, just scan over function and connect points with a line... 
void FiniteFunction::plotFunction(){
  m_function_scan = this->scanFunction(10000);
  m_plotfunction = true;
}

std::vector<double> FiniteFunction::GenerateRandom(int N_random){
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> uniform_dist_x(m_RMin, m_RMax);
  std::uniform_real_distribution<> uniform_dist_T(0, 1);
  double x_i = uniform_dist_x(gen);
  std::vector<double> random_data;
  for (int i = 0; i < N_random; ) {
      std::normal_distribution<> normal_dist(x_i, 0.2);
      double y = normal_dist(gen);
      double A;
      if ((this->callFunction(y)/this->callFunction(x_i))<1){
        A = this->callFunction(y)/this->callFunction(x_i);
      }
      else{
        A=1.;
      }
      double T = uniform_dist_T(gen);
      if (T<A){
        random_data.push_back(y);
        x_i = y;
        i++;
      }
      else{
        continue;
      }
    }
  return random_data;
}

//Transform data points into a format gnuplot can use (histogram) and set flag to enable drawing of data to output plot
//set isdata to true (default) to plot data points in black, set to false to plot sample points in blue
void FiniteFunction::plotData(std::vector<double> &points, int Nbins, bool isdata){
  if (isdata){
    m_data = this->makeHist(points,Nbins);
    m_plotdatapoints = true;
  }
  else{
    m_samples = this->makeHist(points,Nbins);
    m_plotsamplepoints = true;
  }
}


/*
  #######################################################################################################
  ## SUPACPP Note:
  ## The three helper functions below are needed to get the correct format for plotting with gnuplot
  ## In theory you shouldn't have to touch them
  ## However it might be helpful to read through them and understand what they are doing
  #######################################################################################################
 */

//Scan over range of function using range/Nscan steps (just a hack so we can plot the function)
std::vector< std::pair<double,double> > FiniteFunction::scanFunction(int Nscan){
  std::vector< std::pair<double,double> > function_scan;
  double step = (m_RMax - m_RMin)/(double)Nscan;
  double x = m_RMin;
  //We use the integral to normalise the function points
  if (m_Integral == NULL) {
    std::cout << "Integral not set, doing it now" << std::endl;
    this->integral(Nscan);
    std::cout << "integral: " << m_Integral << ", calculated using " << Nscan << " divisions" << std::endl;
  }
  //For each scan point push back the x and y values 
  for (int i = 0; i < Nscan; i++){
    function_scan.push_back( std::make_pair(x,this->callFunction(x)/m_Integral));
    x += step;
  }
  return function_scan;
}

//Function to make histogram out of sampled x-values - use for input data and sampling
std::vector< std::pair<double,double> > FiniteFunction::makeHist(std::vector<double> &points, int Nbins){

  std::vector< std::pair<double,double> > histdata; //Plottable output shape: (midpoint,frequency)
  std::vector<int> bins(Nbins,0); //vector of Nbins ints with default value 0 
  int norm = 0;
  for (double point : points){
    //Get bin index (starting from 0) the point falls into using point value, range, and Nbins
    int bindex = static_cast<int>(floor((point-m_RMin)/((m_RMax-m_RMin)/(double)Nbins)));
    if (bindex<0 || bindex>Nbins){
      continue;
    }
    bins[bindex]++; //weight of 1 for each data point
    norm++; //Total number of data points
  }
  double binwidth = (m_RMax-m_RMin)/(double)Nbins;
  for (int i=0; i<Nbins; i++){
    double midpoint = m_RMin + i*binwidth + binwidth/2; //Just put markers at the midpoint rather than drawing bars
    double normdata = bins[i]/((double)norm*binwidth); //Normalise with N = 1/(Ndata*binwidth)
    histdata.push_back(std::make_pair(midpoint,normdata));
  }
  return histdata;
}

//Function which handles generating the gnuplot output, called in destructor
//If an m_plot... flag is set, the we must have filled the related data vector
//SUPACPP note: They syntax of the plotting code is not part of the course
void FiniteFunction::generatePlot(Gnuplot &gp){

  if (m_plotfunction==true && m_plotdatapoints==true && m_plotsamplepoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 2 lc rgb 'blue' title 'sampled data', '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_samples);
    gp.send1d(m_data);
  }
  else if (m_plotfunction==true && m_plotdatapoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_data);
  }
  else if (m_plotfunction==true && m_plotsamplepoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 2 lc rgb 'blue' title 'sampled data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_samples);
  }
  else if (m_plotfunction==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title 'function'\n";
    gp.send1d(m_function_scan);
  }

  else if (m_plotdatapoints == true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "plot '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_data);
  }

  else if (m_plotsamplepoints == true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "plot '-' with points ps 2 lc rgb 'blue' title 'sampled data'\n";
    gp.send1d(m_samples);
  }
}

std::vector<double> readDataFromFile(const std::string& filePath) {
    std::vector<double> data;
    std::ifstream inputFile(filePath);

    // Check if the file was successfully opened
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open the file: " << filePath << std::endl;
        return data;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        double value;
        if (iss >> value) {
            data.push_back(value);
        } else {
            std::cerr << "Warning: Skipping invalid line: " << line << std::endl;
        }
    }

    inputFile.close();
    return data;
}