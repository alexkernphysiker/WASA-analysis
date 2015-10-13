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
	virtual bool TrackCountTrigger(int CinC,int NinC,int CinF)override;
	virtual bool CentralFirst()override;
	virtual bool ForwardTrackProcessing(WTrack&track)override;
	virtual bool CentralTrackProcessing(WTrack&track)override;
private:
	Analyser2D MissingMass;
	TrackConditionSet Cut,CutFRH1,CutFTH1;
	InterpolationBasedReconstruction He3_theta,He3_phi;
	vector<InterpolationBasedReconstruction> He3_Ekin;
protected:
	TrackDependent Th_m,Th_t,Ph_m,Ph_t,Elow_m,Ehi_m,E_t;
	virtual void debug_cut(WTrack&);
	virtual void debug_notcut(WTrack&);
};
typedef CustomAnalysis<RealData,He3_in_forward> He3Data;
class He3_debug:public He3_in_forward{
public:
	He3_debug();
	virtual ~He3_debug();
protected:
	virtual void debug_cut(WTrack&)override;
	virtual void debug_notcut(WTrack&)override;
private:
	Debug2DSpectraSet Yes,No;
};
class He3eta:public CustomAnalysis<MonteCarlo,He3_debug>{
public:
    He3eta();
};
class He3pi0:public CustomAnalysis<MonteCarlo,He3_debug>{
public:
	He3pi0();
};
#endif 