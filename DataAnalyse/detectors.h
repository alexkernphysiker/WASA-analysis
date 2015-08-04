// this file is distributed under 
// GPL v 3.0 license
#ifndef HGDOYENRKILTMEYD
#define HGDOYENRKILTMEYD
#include <string>
#include <functional>
#include "analysis.h"
#define static_right(A) (static_cast<decltype(A)&&>(A))
enum ForwardDetectorPlane{
	kForwardError=0,
	kFWC1 = 10, kFWC2 = 11, kFTH1 = 1, kFTH2 =2, kFTH3 = 3, 
	kFRH1 = 4, kFRH2 = 5, kFRH3 = 6, kFRH4 = 7,kFRH5 = 8, kFVH =9
};
class ForwardDetectors:public virtual Logger{
public:
	ForwardDetectors();
	virtual ~ForwardDetectors();
protected:
	int ForwadrPlaneCount();
	int ForwardPlaneIndex(ForwardDetectorPlane plane);
	ForwardDetectorPlane ForwadrPlane(int index);
	std::string ForwardPlaneName(int index);
	std::string ForwardPlaneName(ForwardDetectorPlane plane);
	double EDep(WTrack&& track,ForwardDetectorPlane plane);
	double EDep(WTrack&& track,int planeindex);
	double UpperByIndex(int planeindex);
	double Upper(ForwardDetectorPlane plane);
	double ThresholdByIndex(int planeindex);
	double Threshold(ForwardDetectorPlane plane);
	int StopPlaneIndex(WTrack&& track);
	ForwardDetectorPlane StopPlane(WTrack&& track);
	bool UpperThresholdUpTo(WTrack&& track,int planeindex);
private:
	struct plane_data{
	public:
		plane_data(ForwardDetectorPlane p, std::string n,double u,double thr);
		ForwardDetectorPlane plane;
		std::string name;
		double upper;
		double threshold;
	};
	std::vector<plane_data> PlaneData;
};
#endif
