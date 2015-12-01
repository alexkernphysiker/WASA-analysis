// this file is distributed under 
// MIT license
#ifndef mpOBjkfX
#define mpOBjkfX
#include <vector>
#include <string>
#include <functional>
#include <fstream>
#include "math_h/include/exception_math_h.h"
#include "math_h/include/interpolate.h"
#include "config.h"
#include "analysis.h"
#include "log.h"
template<typename... Args>
class InterpolationBasedReconstruction:public virtual Logger{
public:
	typedef std::function<double(Args...)> FUNC;
private:
	std::string m_name;
	bool data_present;
	LinearInterpolation<double> data;
	std::vector<std::pair<double,double>> out;
	FUNC Experiment,Theory;
public:
	InterpolationBasedReconstruction(std::string name,FUNC measured,FUNC theory){
		AddLogSubprefix("InterpolationBasedReconstruction");
		m_name=name;
		AddLogSubprefix(m_name);
		Experiment=measured;
		Theory=theory;
		std::ifstream file;
		file.open((DataFiles+name+".calibration.txt").c_str());
		if(file){
			Log(LogDebug)<<"reading input data";
			double measured,calculated;
			while(file>>measured>>calculated)
				data<<make_pair(measured,calculated);
			file.close();
			data_present=data.size()>0;
		}else{
			data_present=false;
		}
		if(!data_present)Log(NoLog)<<"no input data. Running in simulation mode";
		else Log(NoLog)<<"Input data found. Running in reconstruction mode";
	}
	InterpolationBasedReconstruction(
		const InterpolationBasedReconstruction&source
	){
		m_name=source.m_name;
		data_present=source.data_present;
		data=source.data;
		out=source.out;
		Experiment=source.Experiment;
		Theory=source.Theory;
	}
	virtual ~InterpolationBasedReconstruction(){
		if(!data_present){
			Log(NoLog)<<"saving simulation data";
			std::ofstream file;
			file.open((DataFiles+m_name+".simulation.txt").c_str(),ios_base::app);
			if(file){
				for(auto&p:out)
					file<<p.first<<" "<<p.second<<"\n";
				file.close();
			}
		}
	}
	double Reconstruct(Args... args){
		if(data_present){
			try{
				return data(Experiment(args...));
			}catch(exception){
				Log(LogWarning)<<"Possibly the measured value is out of range.";
				return INFINITY;
			}
		}else{
			out.push_back(make_pair(Experiment(args...),Theory(args...)));
			return INFINITY;
		}
	}
};
#endif 