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
	typedef std::function<double(WTrack&)> delegate;
	InterpolationBasedReconstruction();
	InterpolationBasedReconstruction(
		std::string name,delegate measured,delegate theory
	);
	InterpolationBasedReconstruction(
		const InterpolationBasedReconstruction&source
	);
	virtual ~InterpolationBasedReconstruction();
	double Reconstruct(WTrack&);
private:
	std::string m_name;
	bool data_present;
	LinearInterpolation<double> data;
	std::vector<std::pair<double,double>> out;
	delegate Experiment,Theory;
};
#endif 