// this file is distributed under 
// GPL v 3.0 license
#ifndef hvpdscuo
#define hvpdscuo
#include <string>
#include <functional>
#include <fit.h>
#include <math_h/gnuplot/gnuplot.h>
std::string SimulationDataPath();
void ProcessReconstruction(std::string&&name,double from,double to, std::size_t bins,std::function<bool(double&,double&)> cut,Genetic::RANDOM&engine);
#endif 