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
class InterpolationBasedReconstruction:public virtual Logger{
public:
	typedef std::function<double(WTrack&&)> delegate;
	InterpolationBasedReconstruction(std::string name,delegate measured,delegate theory);
	virtual ~InterpolationBasedReconstruction();
	bool Reconstruct(double&calculated,WTrack&&track);
private:
	std::string m_name;
	bool data_present;
	LinearInterpolation<double> data;
	ofstream out;
	delegate Experiment,Theory;
};
#endif 