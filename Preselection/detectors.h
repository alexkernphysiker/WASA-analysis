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
	int ForwadrPlaneCount()const;
	int ForwardPlaneIndex(ForwardDetectorPlane plane)const;
	ForwardDetectorPlane ForwadrPlane(int index)const;
	std::string&ForwardPlaneName(int index)const;
	std::string&ForwardPlaneName(ForwardDetectorPlane plane)const;
	double EDep(WTrack&track,ForwardDetectorPlane plane)const;
	double EDep(WTrack&track,int planeindex)const;
	double UpperByIndex(int planeindex)const;
	double Upper(ForwardDetectorPlane plane)const;
	double ThresholdByIndex(int planeindex)const;
	double Threshold(ForwardDetectorPlane plane)const;
	int StopPlaneIndex(WTrack&track)const;
	ForwardDetectorPlane StopPlane(WTrack&track)const;
	bool UpperThresholdUpTo(WTrack&track,int planeindex)const;
	void ForwardDetectorTrackMarker(int m,WTrack&track);
private:
	struct plane_data{
	public:
		plane_data(ForwardDetectorPlane p, std::string&&n,double u,double thr);
		ForwardDetectorPlane plane;
		std::string name;
		double upper,threshold;
	};
	std::vector<plane_data> PlaneData;
	std::vector<std::vector<TH2F*>> EDepHist;
};
#endif
