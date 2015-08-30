// this file is distributed under 
// GPL v 3.0 license
#ifndef oVonXYZj
#define oVonXYZj
#include "analysis.h"
#include "detectors.h"
#include "reconstruction.h"
#include "trackprocessing.h"
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
#endif 