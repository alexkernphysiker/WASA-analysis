// this file is distributed under 
// GPL v 3.0 license
#ifndef oVonXYZj
#define oVonXYZj
#include "analysis.h"
#include "detectors.h"
#include "reconstruction.h"
#include "trackprocessing.h"
class He3_at_FRH1:public virtual Analysis,public ForwardDetectors{
public:
	He3_at_FRH1();
	virtual ~He3_at_FRH1();
protected:
	virtual bool EventPreProcessing()override;
	virtual void EventPostProcessing()override;
	virtual bool TrackCountTrigger(int CinC,int NinC,int CinF)override;
	virtual bool CentralFirst()override;
	virtual bool ForwardTrackProcessing(WTrack&&track)override;
	virtual bool CentralTrackProcessing(WTrack&&track)override;
private:
	InterpolationBasedReconstruction He3_Ekin,He3_theta,He3_phi;
	TrackConditionSet Cuts;
	Analyser2D MissingMass;
};
#endif 