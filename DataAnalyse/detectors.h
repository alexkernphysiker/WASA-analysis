#ifndef HGDOYENRKILTMEYD
#define HGDOYENRKILTMEYD
#include <string>
#include <functional>
#include "analysis.h"
enum ForwardDetectorPlane{
	kForwardError=0,
	kFWC1=10,kFWC2=11,kFTH1=1,kFTH2=2,kFTH3=3,
	kFRH1=4,kFRH2=5,kFRH3=6,kFRH4=7,kFRH5=8,kFVH=9
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
	double EDep(WTrack* track,ForwardDetectorPlane plane);
	double EDep(WTrack* track,int planeindex);
	double UpperByIndex(int planeindex);
	double Upper(ForwardDetectorPlane plane);
	double ThresholdByIndex(int planeindex);
	double Threshold(ForwardDetectorPlane plane);
	int StopPlaneIndex(WTrack* track);
	ForwardDetectorPlane StopPlane(WTrack* track);
	bool UpperThresholdUpTo(WTrack* track,int planeindex);
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
template<ParticleType ptype,ForwardDetectorPlane firstplane>
class ForwardDetectorRoutines:public virtual ForwardDetectors{
public:
	typedef std::function<bool(ForwardDetectorPlane,int&)> delegate;
private:
	FDEdep2Ekin *DepKin;
	delegate get_table;
	double coeff;
public:
	ForwardDetectorRoutines(std::string particlename):ForwardDetectors(){
		DepKin = dynamic_cast<FDEdep2Ekin*>(gParameterManager->GetParameterObject("FDEdep2Ekin",particlename.c_str()));
		get_table=[](ForwardDetectorPlane,int&){return true;};
		coeff=1;
	}
	virtual ~ForwardDetectorRoutines(){}
protected:
	void SetGettableFunction(delegate func){
		get_table=func;
	}
	void SetCorrectionCoefficient(double c){
		coeff=c;
	}
	bool ReconstructEkin(WTrack *track,double &Ekin){
		int Edep2Ekin_table=0;
		ForwardDetectorPlane stop_plane=StopPlane(track);
		if(get_table(stop_plane,Edep2Ekin_table)){
			Ekin=DepKin->GetEkin(Edep2Ekin_table,ptype,track->Edep(firstplane,stop_plane),track->Theta());        
			Ekin*=coeff;
			return true;
		}else{
			Ekin=0;
			return false;
		}
	}
};
#endif
