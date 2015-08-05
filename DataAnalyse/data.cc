// this file is distributed under 
// GPL v 3.0 license
#include <fstream>
#include <exception>
#include "data.h"
#include "reconstruction.h"
using namespace std;
RealData::RealData(){
	AddLogSubprefix("Real data analysis");
	fHeader = dynamic_cast<REventHeader*>(gDataManager->GetDataObject("REventHeader","Header"));
	ifstream file;
	file.open((rec_name_prefix+"TIME_IN_CYCLE.txt").c_str());
	if(file.is_open()){
		while(!file.eof()){
			double time,p;
			file>>time>>p;
			p_beam<<make_pair(time/1000.0,p/1000);
		}
		file.close();
	}else{
		Log(LogError)<<"No file with momentum table";
		throw exception();
	}
}
RealData::~RealData(){}
bool RealData::EventProcessingCondition(){
	return true;
}
double RealData::EventWeight(){
	return 1;
}
double RealData::PBeam(){
	double time=fHeader->GetTimeInCycle();
	if((time>=p_beam.min())&&(time<=p_beam.max()))
		return p_beam(time);
	else{
		Log()<<"Time in cycle is out of range where beam energy is defined.";
		return 0;
	}
}
bool RealData::GetTrueParameters(ParticleType type, double& Ekin, double& theta, double& phi){
	return false;
}