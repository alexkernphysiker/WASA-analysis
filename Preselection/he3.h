// this file is distributed under 
// GPL v 3.0 license
#ifndef oVonXYZj
#define oVonXYZj
#include "analysis.h"
#include "detectors.h"
#include "reconstruction.h"
#include "trackprocessing.h"
#include "data.h"
#include "montecarlo.h"
#include "../General/phys_constants.h"
class He3_in_forward:public virtual Analysis,public ForwardDetectors{
public:
	He3_in_forward();
	virtual ~He3_in_forward();
protected:
	virtual bool EventPreProcessing()override;
	virtual void EventPostProcessing()override;
	virtual bool TrackCountTrigger(int CinC,int NinC,int CinF,int NinF)override;
private:
	TrackConditionSet Cut,CutFTH1,CutFRH1;
	Analyser2D MissingMass;
	InterpolationBasedReconstruction He3_theta,He3_phi;
	vector<InterpolationBasedReconstruction> He3_Ekin;
	TrackDependent Th_m,Th_t,Ph_m,Ph_t,Elow_m,Ehi_m,E_t;
};
typedef CustomAnalysis<RealData,He3_in_forward> He3Data;
class He3eta:public CustomAnalysis<MonteCarlo,He3_in_forward>{
public:
	He3eta();
};
class He3pi0:public CustomAnalysis<MonteCarlo,He3_in_forward>{
public:
	He3pi0();
};
#endif 