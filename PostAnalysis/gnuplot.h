#ifndef VIJVUSSC
#define VIJVUSSC
#include <vector>
#include <functional>
#include <memory>
#include <fitpoints.h>
#include <math_h/interpolate.h>
#include "gethist.h"
void SetPlotOutput(std::string outpath);
class Plot{
public:
	Plot();
	~Plot();
	Plot &operator<<(std::string line);
	Plot &Points(std::string name,std::shared_ptr<Fit::FitPoints> points);
	Plot &Points(std::string name,LinearInterpolation<double> &points);
	Plot &Points(std::string name,LinearInterpolation<double> &points,std::function<double(double)> error);
	Plot &Function(std::string name,std::function<double(double)> func,double from,double to,double step);
	Plot &Function(std::string name,std::function<double(double)> func,std::function<double(double)> error,double from,double to,double step);
private:
	std::vector<std::string> lines;
	std::vector<std::string> plots;
};
#endif