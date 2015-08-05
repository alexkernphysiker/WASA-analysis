// this file is distributed under 
// GPL v 3.0 license
#include <fstream>
#include "reconstruction.h"
using namespace std;
InterpolationBasedReconstruction::InterpolationBasedReconstruction(
	std::string name,delegate measured,delegate theory,
	double from,double to, int bins
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
	if(!data_present){
		Log(NoLog)<<"no input data. Running in simulation mode";
		output=new TH2F(name.c_str(),"",bins,from,to,bins,from,to);
		gHistoManager->Add(output,"Reconstruction");
	}else Log(NoLog)<<"Input data found. Running in reconstruction mode";
}
InterpolationBasedReconstruction::~InterpolationBasedReconstruction(){}
bool InterpolationBasedReconstruction::Reconstruct(double& calculated,WTrack&&track){
	if(data_present){
		try{
			calculated=data(Experiment(static_cast<WTrack&&>(track)));
			return true;
		}catch(exception){
			Log(LogWarning)<<"Possibly the measured value is out of range.";
			return false;
		}
	}else{
		output->Fill(Experiment(static_cast<WTrack&&>(track)),Theory(static_cast<WTrack&&>(track)));
		return false;
	}
}

