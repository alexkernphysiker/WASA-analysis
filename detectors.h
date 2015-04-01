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
template<ParticleType ptype>
class ForwardDetectorRoutines:public virtual ForwardDetectors{
private:
	FDEdep2Ekin *DepKin;
	std::function<bool(int,int&,int&)> check;
	double coeff;
public:
	ForwardDetectorRoutines(string tablename,string particlename,double corr_coef):ForwardDetectors(){
		DepKin = dynamic_cast<FDEdep2Ekin*>(gParameterManager->GetParameterObject(tablename.c_str(),particlename.c_str()));
		check=[](int,int&,int&){return true;};
		coeff=corr_coef;
	}
	virtual ~ForwardDetectorRoutines(){}
protected:
	void SetReconstructionCheckFunction(std::function<bool(int,int&,int&)> func){
		check=func;
	}
	ForwardDetectorPlane StoppingPlane(WTrack *track){
		return DepKin->GetStoppingPlane(track);
	}
	bool ReconstructEkin(WTrack *track,double &Ekin){
		Int_t Edep2Ekin_table=0;
		Int_t last_plane=0;
		Int_t StopPlane=StoppingPlane(track);
		if(check(StopPlane,last_plane,Edep2Ekin_table)){
			Ekin=DepKin->GetEkin(Edep2Ekin_table,14,track->Edep(4,last_plane),track->Theta()); //kinetic energy of 3He	        
			Ekin*=coeff;
			return true;
		}else
			return false;
	}
};
#endif
