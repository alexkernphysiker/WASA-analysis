// this file is distributed under 
// GPL v 3.0 license
#ifndef oVonXYZj
#define oVonXYZj
#include "analysis.h"
#include "detectors.h"
#include "reconstruction.h"
#include "trackprocessing.h"
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
	TrackConditionSet Cut,CutFRH1,CutFTH1;
	InterpolationBasedReconstruction He3_Ekin,He3_theta,He3_phi;
	Analyser2D MissingMass;
};
class He3eta:public He3_in_forward{
public:
    He3eta();
};
template<unsigned int count>
class He3pi0:public He3_in_forward{
public:
	He3pi0():He3_in_forward(){
		for(unsigned int i=0;i<count;i++)
			AddParticleToFirstVertex(kPi0,m_pi0);
	}
};
#endif 