// this file is distributed under 
// GPL v 3.0 license
#include <vector>
#include <utility>
#include "detectors.h"
using namespace std;
ForwardDetectors::plane_data::plane_data(ForwardDetectorPlane p, string n, double u,double thr)
	:plane(p),name(n),upper(u),threshold(thr){}
ForwardDetectors::ForwardDetectors(){
	AddLogSubprefix("Forward detector");
	PlaneData.push_back(plane_data(kFWC1,"FWC1",0.03,0.002));
	PlaneData.push_back(plane_data(kFWC2,"FWC2",0.03,0.002));
	//PlaneData.push_back(plane_data(kFPC,"FPC",0.03,0.002));
	PlaneData.push_back(plane_data(kFTH1,"FTH1",0.05,0.0015));
	PlaneData.push_back(plane_data(kFRH1,"FRH1",0.3 ,0.001));
	PlaneData.push_back(plane_data(kFRH2,"FRH2",0.3 ,0.001));
}
ForwardDetectors::~ForwardDetectors(){}
int ForwardDetectors::ForwadrPlaneCount(){
	return PlaneData.size();
}
int ForwardDetectors::ForwardPlaneIndex(ForwardDetectorPlane plane){
	for(size_t i=0;i<PlaneData.size();i++)
		if(PlaneData[i].plane==plane)
			return i;
		Log()<<"forward plane not found";
	return -1;
}
ForwardDetectorPlane ForwardDetectors::ForwadrPlane(int index){
	if((index>=0)&&(index<int(PlaneData.size())))
		return PlaneData[index].plane;
	Log()<<"forward plane by index not found";
	return kForwardError;
}
string ForwardDetectors::ForwardPlaneName(int index){
	if((index>=0)&&(index<int(PlaneData.size())))
		return PlaneData[index].name;
	Log()<<"forward plane by index not found";
	return "<NoPlane>";
}
string ForwardDetectors::ForwardPlaneName(ForwardDetectorPlane plane){
	return ForwardPlaneName(ForwardPlaneIndex(plane));
}
double ForwardDetectors::EDep(WTrack&& track, ForwardDetectorPlane plane){
	return track.Edep(plane);
}

double ForwardDetectors::EDep(WTrack&& track, int planeindex){
	return EDep(static_right(track),ForwadrPlane(planeindex));
}
double ForwardDetectors::UpperByIndex(int planeindex){
	return PlaneData[planeindex].upper;
}
double ForwardDetectors::Upper(ForwardDetectorPlane plane){
	return PlaneData[ForwardPlaneIndex(plane)].upper;
}
double ForwardDetectors::ThresholdByIndex(int planeindex){
	return PlaneData[planeindex].threshold;
}
double ForwardDetectors::Threshold(ForwardDetectorPlane plane){
	return PlaneData[ForwardPlaneIndex(plane)].threshold;
}
bool ForwardDetectors::UpperThresholdUpTo(WTrack&& track,int planeindex){
	SubLog log=Log(LogDebug);
	log<<"UpperThresholdUpTo call";
	bool res=true;
	for(int c=0;c<=planeindex;c++)
		res&=(EDep(static_right(track),c)>ThresholdByIndex(c));
	if(res)log<<"YES";else log<<"NO";
	return res;
}
int ForwardDetectors::StopPlaneIndex(WTrack&& track){
	int res=-1;
	for(size_t i=0;i<PlaneData.size();i++)
		if(EDep(static_right(track),i)>PlaneData[i].threshold){
			if(res==(int(i)-1))res++;
			else res=-1;
		}
	return res;
}
ForwardDetectorPlane ForwardDetectors::StopPlane(WTrack&& track){
	return ForwadrPlane(StopPlaneIndex(static_right(track)));
}
