// this file is distributed under 
// MIT license
#ifndef RLDdVIXV
#define RLDdVIXV
#include <memory>
#include <string>
#include <functional>
#include <fit.h>
std::shared_ptr<Genetic::FitPoints> ReadWeightedFrom2D(
	std::string name,
	double fromX, double toX, unsigned int binsX,
	double fromY, double toY, unsigned int binsY,
	std::function<bool(double&,double&)> processing
);
#endif 