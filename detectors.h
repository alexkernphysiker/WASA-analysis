#ifndef HGDOYENRKILTMEYD
#define HGDOYENRKILTMEYD
#include <string>
#include "analysis.h"
enum ForwardDetectorPlane{
	kForwardError=0,
	kFWC1=10,kFWC2=11,kFTH1=1,kFTH2=2,kFTH3=3,
	kFRH1=4,kFRH2=5,kFRH3=6,kFRH4=7,kFRH5=8,kFVH=9
};
class ForwardDetectors{
public:
	ForwardDetectors();
	virtual ~ForwardDetectors();
protected:
	int ForwadrPlaneCount();
	int ForwardPlaneIndex(ForwardDetectorPlane plane);
	ForwardDetectorPlane ForwadrPlane(int index);
	std::string ForwardPlaneName(int index);
	std::string ForwardPlaneName(ForwardDetectorPlane plane);
	double EDep(WTrack* track,int planeindex);
	bool IsForwardPlaneDep(WTrack* track,int planeindex);
	bool IsForwardPlaneDep(WTrack* track,ForwardDetectorPlane plane);
	int ForwardStopPlaneIndex(WTrack* track);
	ForwardDetectorPlane ForwardStopPlane(WTrack* track);
private:
	struct plane_data{
	public:
		plane_data(ForwardDetectorPlane p,double t, std::string n);
		ForwardDetectorPlane plane;
		double threshold;
		std::string name;
	};
	std::vector<plane_data> PlaneData;
};
#endif
