// this file is distributed under 
// MIT license
#ifndef oVonXYZj
#define oVonXYZj
#include "analysis.h"
#include "detectors.h"
#include "reconstruction.h"
#include "trackprocessing.h"
#include "data.h"
#include "montecarlo.h"
#include "../phys_constants.h"
class He3_in_forward:public virtual Analysis,public ForwardDetectors{
public:
	He3_in_forward(double Q_lo,double Q_hi,unsigned int bins);//Q is in MeV
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
class He3_Modification_for_reconstruction:public He3_in_forward{
public:
	He3_Modification_for_reconstruction();
};
class He3_Modification_for_eta:public He3_in_forward{
public:
	He3_Modification_for_eta();
};
template<class He3>class He3eta:public CustomAnalysis<MonteCarlo,He3>{
	public:He3eta(){He3_in_forward::AddParticleToFirstVertex(kEta,m_eta);}
};
template<class He3>class He3pi0:public CustomAnalysis<MonteCarlo,He3>{
	public:He3pi0(){He3_in_forward::AddParticleToFirstVertex(kPi0,m_pi0);}
};


typedef CustomAnalysis<MonteCarlo,He3eta<He3_Modification_for_reconstruction>> RE_He3eta;
typedef CustomAnalysis<MonteCarlo,He3pi0<He3_Modification_for_reconstruction>> RE_He3pi0;
typedef CustomAnalysis<MonteCarlo,He3eta<He3_Modification_for_eta>> MC_He3eta;
typedef CustomAnalysis<MonteCarlo,He3pi0<He3_Modification_for_eta>> MC_He3pi0;
typedef CustomAnalysis<RealData,He3_Modification_for_eta> Data_He3;
#endif 