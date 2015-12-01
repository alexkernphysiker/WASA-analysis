// this file is distributed under 
// MIT license
#include <vector>
#include <utility>
#include "detectors.h"
using namespace std;
ForwardDetectors::plane_data::plane_data(ForwardDetectorPlane p, string&&n, double u,double thr)
	:plane(p),name(n),upper(u),threshold(thr){}
ForwardDetectors::ForwardDetectors(int markers){
	AddLogSubprefix("Forward detector");
	PlaneData.push_back(plane_data(kFWC1,"FWC1",0.03,0.002));
	PlaneData.push_back(plane_data(kFWC2,"FWC2",0.03,0.002));
	//PlaneData.push_back(plane_data(kFPC,"FPC",0.03,0.002));
	PlaneData.push_back(plane_data(kFTH1,"FTH1",0.05,0.0015));
	PlaneData.push_back(plane_data(kFRH1,"FRH1",0.5 ,0.001));
	PlaneData.push_back(plane_data(kFRH2,"FRH2",0.5 ,0.001));
	
	for(int i=0,n=ForwadrPlaneCount()-1;i<n;i++){
		string histname=ForwardPlaneName(i)+"_vs_"+ForwardPlaneName(i+1)+"_";
		EDepHist.push_back({});
		for(int m=0;m<markers;m++){
			string thishistname=histname+to_string(m);
			TH2F* hist=new TH2F(thishistname.c_str(),"",128,0,UpperByIndex(i+1),128,0,UpperByIndex(i));
			EDepHist[i].push_back(hist);
			gHistoManager->Add(hist,"ForwardDetector_EDep");
		}
	}
	
}
ForwardDetectors::~ForwardDetectors(){}

int ForwardDetectors::ForwadrPlaneCount()const{return PlaneData.size();}
int ForwardDetectors::ForwardPlaneIndex(ForwardDetectorPlane plane)const{
	for(size_t i=0;i<PlaneData.size();i++)
		if(PlaneData[i].plane==plane)
			return i;
	return -1;
}
ForwardDetectorPlane ForwardDetectors::ForwadrPlane(int index)const{
	if((index>=0)&&(index<int(PlaneData.size())))
		return PlaneData[index].plane;
	return kForwardError;
}
string&ForwardDetectors::ForwardPlaneName(int index)const{
	if((index>=0)&&(index<int(PlaneData.size())))
		return const_cast<string&>(PlaneData[index].name);
	static string noplane="<NoPlane>";
	return const_cast<string&>(noplane);
}
string&ForwardDetectors::ForwardPlaneName(ForwardDetectorPlane plane)const{
	return ForwardPlaneName(ForwardPlaneIndex(plane));
}
double ForwardDetectors::EDep(WTrack&track, ForwardDetectorPlane plane)const{
	return track.Edep(plane);
}
double ForwardDetectors::EDep(WTrack&track, int planeindex)const{
	return EDep(track,ForwadrPlane(planeindex));
}
double ForwardDetectors::UpperByIndex(int planeindex)const{
	return PlaneData[planeindex].upper;
}
double ForwardDetectors::Upper(ForwardDetectorPlane plane)const{
	return PlaneData[ForwardPlaneIndex(plane)].upper;
}
double ForwardDetectors::ThresholdByIndex(int planeindex)const{
	return PlaneData[planeindex].threshold;
}
double ForwardDetectors::Threshold(ForwardDetectorPlane plane)const{
	return PlaneData[ForwardPlaneIndex(plane)].threshold;
}
bool ForwardDetectors::UpperThresholdUpTo(WTrack&track,int planeindex)const{
	bool res=true;
	for(int c=0;c<=planeindex;c++)
		res&=(EDep(track,c)>ThresholdByIndex(c));
	return res;
}
int ForwardDetectors::StopPlaneIndex(WTrack&track)const{
	int res=-1;
	for(size_t i=0;i<PlaneData.size();i++)
		if(EDep(track,i)>PlaneData[i].threshold){
			if(res==(int(i)-1))res++;
			else res=-1;
		}
	return res;
}
ForwardDetectorPlane ForwardDetectors::StopPlane(WTrack&track)const{
	return ForwadrPlane(StopPlaneIndex(track));
}

void ForwardDetectors::ForwardDetectorTrackMarker(int m,WTrack&track){
	for(int i=0,n=ForwadrPlaneCount()-1;i<n;i++)
		EDepHist[i][m]->Fill(EDep(track,i+1),EDep(track,i));
}
