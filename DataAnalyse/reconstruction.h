// this file is distributed under 
// GPL v 3.0 license
#ifndef mpOBjkfX
#define mpOBjkfX
#include <vector>
#include <string>
#include <functional>
#include "math_h/interpolate.h"
#include "analysis.h"
#include "log.h"
const string rec_name_prefix="../Reconstruction/";
class InterpolationBasedReconstruction:public virtual Logger{
public:
	typedef std::function<double(WTrack&&)> delegate;
	InterpolationBasedReconstruction(
		std::string name,delegate measured,delegate theory,
		double from,double to, int bins
	);
	virtual ~InterpolationBasedReconstruction();
	bool Reconstruct(double&calculated,WTrack&&track);
private:
	std::string m_name;
	bool data_present;
	LinearInterpolation<double> data;
	TH2F* output;
	delegate Experiment,Theory;
};
#endif 