#include <vector>
#include <utility>
#include "detectors.h"
using namespace std;
ForwardDetectors::plane_data::plane_data(ForwardDetectorPlane p, string n, double u,double thr)
	:plane(p),name(n),upper(u),threshold(thr){}
ForwardDetectors::ForwardDetectors(){
	PlaneData.push_back(plane_data(kFWC1,"FWC1",0.03,0.0005));
	PlaneData.push_back(plane_data(kFWC2,"FWC2",0.03,0.0005));
	PlaneData.push_back(plane_data(kFTH1,"FTH1",0.05,0.0005));
	PlaneData.push_back(plane_data(kFTH2,"FTH2",0.05,0.0005));
	PlaneData.push_back(plane_data(kFTH3,"FTH3",0.05,0.0005));
	PlaneData.push_back(plane_data(kFRH1,"FRH1",0.3 ,0.0005));
	PlaneData.push_back(plane_data(kFRH2,"FRH2",0.3 ,0.0005));
	PlaneData.push_back(plane_data(kFRH3,"FRH3",0.3 ,0.0005));
	PlaneData.push_back(plane_data(kFRH4,"FRH4",0.3 ,0.0005));
	PlaneData.push_back(plane_data(kFRH5,"FRH5",0.3 ,0.0005));
	PlaneData.push_back(plane_data(kFVH ,"FVH" ,0.05,0.0005));
}
ForwardDetectors::~ForwardDetectors(){}
int ForwardDetectors::ForwadrPlaneCount(){
	return PlaneData.size();
}
int ForwardDetectors::ForwardPlaneIndex(ForwardDetectorPlane plane){
	for(size_t i=0;i<PlaneData.size();i++)
		if(PlaneData[i].plane==plane)
			return i;
	throw;//no such detector here
}
ForwardDetectorPlane ForwardDetectors::ForwadrPlane(int index){
	return PlaneData[index].plane;
}
string ForwardDetectors::ForwardPlaneName(int index){
	return PlaneData[index].name;
}
string ForwardDetectors::ForwardPlaneName(ForwardDetectorPlane plane){
	return PlaneData[ForwardPlaneIndex(plane)].name;
}
double ForwardDetectors::EDep(WTrack* track, int planeindex){
	return track->Edep(PlaneData[planeindex].plane);
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
bool ForwardDetectors::ThresholdCondition(WTrack* track,ForwardDetectorPlane plane){
	int index=ForwardPlaneIndex(plane);
	bool res=true;
	for(int c=0;c<index;c++)
		res&=(EDep(track,c)>ThresholdByIndex(c));
	return res;
}

