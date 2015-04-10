#include <fstream>
#include <iostream>
#include "data.h"
using namespace std;
RealData::RealData(){
	fHeader = dynamic_cast<REventHeader*>(gDataManager->GetDataObject("REventHeader","Header"));
	ifstream file;
	file.open("TIME_IN_CYCLE_MAY2014.dat");
	if(file.is_open()){
		while(!file.eof()){
			double time,p;
			file>>time>>p;
			p_beam<<make_pair(time/1000.0,p);
		}
		file.close();
	}else{
		printf(">>>>>>>No beam momenta table file\n");
		throw;
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
	return p_beam(fHeader->GetTimeInCycle());
}
void RealData::PrepareCheck(){}
void RealData::CheckParticleTrack(ParticleType type, double Ekin, double theta, double phi){}
