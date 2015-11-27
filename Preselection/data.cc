// this file is distributed under 
// GPL v 3.0 license
#include <fstream>
#include <exception>
#include "config.h"
#include "data.h"
#include "reconstruction.h"
using namespace std;
RealData::RealData():BeamMomenta("Time.2.PBeam",
		[this](){return fHeader->GetTimeInCycle()*1000.0;},
		[](){return INFINITY;}
){
	AddLogSubprefix("Real data analysis");
	fHeader = dynamic_cast<REventHeader*>(gDataManager->GetDataObject("REventHeader","Header"));
}
RealData::~RealData(){}
bool RealData::EventProcessingCondition(){
	return true;
}
void RealData::PrepairForEventAnalysis(){
	CachePBeam(BeamMomenta.Reconstruct()/1000.0);
}