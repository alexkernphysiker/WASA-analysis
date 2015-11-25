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
#include "trackprocessing.h"
class InterpolationBasedReconstruction:public virtual Logger{
public:
	InterpolationBasedReconstruction();
	InterpolationBasedReconstruction(
		std::string name,ValueTrackDependent measured,ValueTrackDependent theory
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
	ValueTrackDependent Experiment,Theory;
};
#endif 