// this file is distributed under 
// GPL v 3.0 license
#include <fstream>
#include "reconstruction.h"
using namespace std;
const string nameprefix="../Reconstruction";
InterpolationBasedReconstruction::InterpolationBasedReconstruction(string name,delegate measured,delegate theory){
	AddLogSubprefix("InterpolationBasedReconstruction");
	m_name=name;
	AddLogSubprefix("name");
	Experiment=measured;
	Theory=theory;
	ifstream file;
	file.open((nameprefix+name+".calibration.txt").c_str());
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
	if(!data_present)Log(LogDebug)<<"no input data. Simulation mode";
}
InterpolationBasedReconstruction::~InterpolationBasedReconstruction(){
	if(!data_present){
		Log(LogDebug)<<"saving simulation data";
		ofstream file;
		file.open((nameprefix+m_name+".simulation.txt").c_str());
		if(file.is_open()){
			for(auto&p:data)
				file<<p.first<<" "<<p.second<<"\n";
			file.close();
		}else{
			Log(LogError)<<"cannot write output file.";
		}
	}
}
bool InterpolationBasedReconstruction::Reconstruct(double& calculated,WTrack&&track){
	if(data_present){
		try{
			calculated=data(Experiment(static_cast<WTrack&&>(track)));
			return true;
		}catch(exception){
			Log(LogWarning)<<"error during reconstruction. Possibly the measured value is out of range.";
			return false;
		}
	}else{
		data<<make_pair(Experiment(static_cast<WTrack&&>(track)),Theory(static_cast<WTrack&&>(track)));
		return false;
	}
}

