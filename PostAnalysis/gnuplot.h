#ifndef VIJVUSSC
#define VIJVUSSC
#include <vector>
#include <functional>
#include <memory>
#include <fitpoints.h>
#include "gethist.h"
class Plot{
public:
	Plot(std::string out);
	virtual ~Plot();
	Plot &Points(std::string file,std::shared_ptr<Fit::FitPointsAbstract> points);
	Plot &Hist(std::string file,hist &points);
	Plot &Function(std::string file,std::function<double(double)> func,double from,double to,double step);
	void Out(std::string scriptname,bool show=false);
private:
	std::string outpath;
	std::vector<std::string> lines;
};
#endif