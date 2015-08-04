// this file is distributed under 
// GPL v 3.0 license
#ifndef EVHZDETEOYGNMFFP
#define EVHZDETEOYGNMFFP
#include "analysis.h"
#include "detectors.h"
#include "reconstruction.h"
class He3eta_gg_:public virtual Analysis,public ForwardDetectors{
public:
	He3eta_gg_();
	virtual ~He3eta_gg_();
protected:
	virtual bool EventPreProcessing(TVector3 &&pbeam)override;
	virtual void EventPostProcessing(TVector3 &&pbeam)override;
	virtual bool TrackCountTrigger(int CinC,int NinC,int CinF)override;
	virtual bool CentralFirst()override;
	virtual bool ForwardTrackProcessing(WTrack&&track,TVector3&&pbeam)override;
	virtual bool CentralTrackProcessing(WTrack&&track,TVector3&&pbeam)override;
private:
	InterpolationBasedReconstruction He3_Ekin,He3_theta,He3_phi;
	vector<TH2F*> EDepHist;
	vector<TH2F*> EDepFilteredHist;
	vector<TH1F*> MissingMassDetailed;
	TH1F *MissingMass,*DependenceOnPBeam,*P_Beam;
};
#endif