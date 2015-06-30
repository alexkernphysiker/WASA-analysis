#ifndef EVHZDETEOYGNMFFP
#define EVHZDETEOYGNMFFP
#include "analysis.h"
#include "detectors.h"
class He3eta_gg_:public virtual Analysis
	// Though the table used here contains data about He3 it has
	// the descriptor number like for protons. 
	// It's a mistake made in the table, not in my program
	,public ForwardDetectorRoutines<kProton,kFRH1>
{
public:
	He3eta_gg_();
	virtual ~He3eta_gg_();
protected:
	virtual bool EventPreProcessing(TVector3 &&pbeam)override;
	virtual void EventPostProcessing(TVector3 &&pbeam)override;
	virtual bool TrackCountTrigger(int CinC,int NinC,int CinF)override;
	virtual bool CentralFirst()override;
	virtual bool ForwardTrackProcessing(WTrack&& track,TVector3 &&pbeam)override;
	virtual bool CentralTrackProcessing(WTrack&& track,TVector3 &&pbeam)override;
private:
	vector<TH2F*> EDepHist;
	vector<TH2F*> EDepFilteredHist;
	vector<TH1F*> MissingMassDetailed;
	TH1F *MissingMass,*DependenceOnPBeam,*P_Beam;
};
#endif