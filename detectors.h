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
	double Upper(int planeindex);
	double Upper(ForwardDetectorPlane plane);
private:
	struct plane_data{
	public:
		plane_data(ForwardDetectorPlane p,double u, std::string n);
		ForwardDetectorPlane plane;
		double upper;
		std::string name;
	};
	std::vector<plane_data> PlaneData;
};
template<ParticleType ptype,ForwardDetectorPlane firstplane>
class ForwardDetectorRoutines:public virtual ForwardDetectors{
private:
	FDEdep2Ekin *DepKin;
	std::function<bool(int&,int&)> get_table;
	double coeff;
public:
	ForwardDetectorRoutines(std::string particlename):ForwardDetectors(){
		DepKin = dynamic_cast<FDEdep2Ekin*>(gParameterManager->GetParameterObject("FDEdep2Ekin",particlename.c_str()));
		g=[](int&,int&){return true;};
		coeff=1;
	}
	virtual ~ForwardDetectorRoutines(){}
protected:
	void SetReconstructionGettableFunction(std::function<bool(int,int&,int&)> func){
		get_table=func;
	}
	void SetCorrectionCoefficient(double c){
		coeff=c;
	}
	int StoppingPlane(WTrack *track){
		return DepKin->GetStoppingPlane(track);
	}
	bool ReconstructEkin(WTrack *track,double &Ekin){
		int Edep2Ekin_table=0;
		int stop_plane=StoppingPlane(track);
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
