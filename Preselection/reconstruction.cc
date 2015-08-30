// this file is distributed under 
// GPL v 3.0 license
#include <fstream>
#include "reconstruction.h"
using namespace std;
InterpolationBasedReconstruction::InterpolationBasedReconstruction(
	std::string name,delegate measured,delegate theory
){
	AddLogSubprefix("InterpolationBasedReconstruction");
	m_name=name;
	AddLogSubprefix(m_name);
	Experiment=measured;
	Theory=theory;
	ifstream file;
	file.open((rec_name_prefix+name+".calibration.txt").c_str());
	if(file.is_open()){
		Log(LogDebug)<<"reading input data";
		while(!file.eof()){
			double measured,calculated;
			file>>measured>>calculated;
			data<<make_pair(measured,calculated);
		}
		file.close();
		data_present=data.size()>0;
	}else{
		data_present=false;
	}
	if(!data_present)Log(NoLog)<<"no input data. Running in simulation mode";
	else Log(NoLog)<<"Input data found. Running in reconstruction mode";
}
InterpolationBasedReconstruction::~InterpolationBasedReconstruction(){
	if(!data_present){
		Log(NoLog)<<"saving simulation data";
		ofstream file;
		file.open((rec_name_prefix+m_name+".simulation.txt").c_str());
		if(file.is_open()){
			for(auto&p:out)
				file<<p.first<<" "<<p.second<<"\n";
			file.close();
		}
	}
}
double InterpolationBasedReconstruction::Reconstruct(WTrack&track){
	if(data_present){
		try{
			return data(Experiment(track));
		}catch(exception){
			Log(LogWarning)<<"Possibly the measured value is out of range.";
			return INFINITY;
		}
	}else{
		out.push_back(make_pair(Experiment(track),Theory(track)));
		return INFINITY;
	}
}

