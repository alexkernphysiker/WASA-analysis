#include <vector>
#include <utility>
#include "detectors.h"
using namespace std;
ForwardDetectors::plane_data::plane_data(ForwardDetectorPlane p, double t, double u, string n)
	:plane(p),threshold(t),upper(u),name(n){}
ForwardDetectors::ForwardDetectors(){
	PlaneData.push_back(plane_data(kFWC1,0.0025,0.25,"FWC1"));
	PlaneData.push_back(plane_data(kFWC2,0.0025,0.05 ,"FWC2"));
	PlaneData.push_back(plane_data(kFTH1,0.0025,0.05,"FTH1"));
	PlaneData.push_back(plane_data(kFTH2,0.0025,0.05,"FTH2"));
	PlaneData.push_back(plane_data(kFTH3,0.0025,0.05,"FTH3"));
	PlaneData.push_back(plane_data(kFRH1,0.004 ,0.3 ,"FRH1"));
	PlaneData.push_back(plane_data(kFRH2,0.0025,0.3 ,"FRH2"));
	PlaneData.push_back(plane_data(kFRH3,0.0025,0.3 ,"FRH3"));
	PlaneData.push_back(plane_data(kFRH4,0.0035,0.3 ,"FRH4"));
	PlaneData.push_back(plane_data(kFRH5,0.004 ,0.3 ,"FRH5"));
	PlaneData.push_back(plane_data(kFVH ,0.004 ,0.05,"FVH" ));
}
ForwardDetectors::~ForwardDetectors(){}
int ForwardDetectors::ForwadrPlaneCount(){
	return PlaneData.size();
}
int ForwardDetectors::ForwardPlaneIndex(ForwardDetectorPlane plane){
	for(int i=0;i<PlaneData.size();i++)
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
double ForwardDetectors::Upper(int planeindex){
	return PlaneData[planeindex].upper;
}
double ForwardDetectors::Upper(ForwardDetectorPlane plane){
	return PlaneData[ForwardPlaneIndex(plane)].upper;
}
bool ForwardDetectors::IsForwardPlaneDep(WTrack* track,int planeindex){
	return track->Edep(PlaneData[planeindex].plane)>PlaneData[planeindex].threshold;
}
bool ForwardDetectors::IsForwardPlaneDep(WTrack* track,ForwardDetectorPlane plane){
	return IsForwardPlaneDep(track,ForwardPlaneIndex(plane));
}
int ForwardDetectors::ForwardStopPlaneIndex(WTrack* track){
	int res=-1;int cnt=-1;
	for(auto planedata:PlaneData){
		cnt++;
		if(track->Edep(planedata.plane)>planedata.threshold){
			if(res<(cnt-1))
				return -1;
			res++;
		}
	}
	return res;
}
ForwardDetectorPlane ForwardDetectors::ForwardStopPlane(WTrack* track){
	int res=-1;int cnt=-1;
	for(auto planedata:PlaneData){
		cnt++;
		if(track->Edep(planedata.plane)>planedata.threshold){
			if(res<(cnt-1))
				return kForwardError;
			res++;
		}
	}
	return PlaneData[res].plane;
}

