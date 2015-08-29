// this file is distributed under 
// GPL v 3.0 license
#ifndef HGDOYENRKILTMEYD
#define HGDOYENRKILTMEYD
#include <string>
#include <functional>
#include <TH2F.h>
#include "analysis.h"
enum ForwardDetectorPlane{
	kForwardError=0,
	kFWC1 = 10, kFWC2 = 11, kFTH1 = 1, kFTH2 =2, kFTH3 = 3, 
	kFRH1 = 4, kFRH2 = 5, kFRH3 = 6, kFRH4 = 7,kFRH5 = 8, kFVH =9
};
class ForwardDetectors:public virtual Logger{
public:
	ForwardDetectors(int markers);
	virtual ~ForwardDetectors();
protected:
	int ForwadrPlaneCount();
	int ForwardPlaneIndex(ForwardDetectorPlane plane);
	ForwardDetectorPlane ForwadrPlane(int index);
	std::string&&ForwardPlaneName(int index);
	std::string&&ForwardPlaneName(ForwardDetectorPlane plane);
	double EDep(WTrack&track,ForwardDetectorPlane plane);
	double EDep(const WTrack&track,int planeindex);
	double UpperByIndex(int planeindex);
	double Upper(ForwardDetectorPlane plane);
	double ThresholdByIndex(int planeindex);
	double Threshold(ForwardDetectorPlane plane);
	int StopPlaneIndex(const WTrack&track);
	ForwardDetectorPlane StopPlane(const WTrack&track);
	bool UpperThresholdUpTo(const WTrack&track,int planeindex);
	void ForwardDetectorTrackMarker(int m,const WTrack&track);
private:
	struct plane_data{
	public:
		plane_data(ForwardDetectorPlane p, std::string n,double u,double thr);
		ForwardDetectorPlane plane;
		std::string name;
		double upper,threshold;
	};
	std::vector<plane_data> PlaneData;
	std::vector<std::vector<TH2F*>> EDepHist;
};
#endif
