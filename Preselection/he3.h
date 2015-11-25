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
class He3_in_forward:public virtual Analysis,public ForwardDetectors{
public:
	He3_in_forward();
	virtual ~He3_in_forward();
protected:
	virtual bool EventPreProcessing()override;
	virtual void EventPostProcessing()override;
	virtual bool TrackCountTrigger(int CinC,int NinC,int CinF,int NinF)override;
	void AddCondition(ConditionTrackDependent condition);
private:
	vector<TrackConditionSet> ForwardLayerCuts;
	TrackConditionSet Reconstruction;
	Analyser2D MissingMass;
	vector<ConditionTrackDependent> AdditionalConditions;
};
//Here will be different modifications of analysis
class He3_Modification_for_eta:public He3_in_forward{public:He3_Modification_for_eta();};

typedef He3_in_forward He3Catcher;//This analysis class will be used for all procedures
typedef CustomAnalysis<RealData,He3Catcher> Data_He3;
class MC_He3eta:public CustomAnalysis<MonteCarlo,He3Catcher>{public:MC_He3eta();};
class MC_He3pi0:public CustomAnalysis<MonteCarlo,He3Catcher>{public:MC_He3pi0();};
#endif 