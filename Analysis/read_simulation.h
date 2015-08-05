// this file is distributed under 
// GPL v 3.0 license
#ifndef RLDdVIXV
#define RLDdVIXV
#include <memory>
#include <string>
#include <fit.h>
std::shared_ptr<Genetic::FitPoints> ReadWeightedFrom2D(std::string name,double from, double to, unsigned int bins);
#endif 