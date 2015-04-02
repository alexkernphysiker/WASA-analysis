#ifndef EVHZDETEOYGNMFFP
#define EVHZDETEOYGNMFFP
#include "analysis.h"
#include "detectors.h"
class He3eta:public virtual Analysis
	//Protony bo ta tabelka jest trochę dziwna
	,public ForwardDetectorRoutines<kProton,kFRH1>
{
public:
	He3eta();
	virtual ~He3eta();
protected:
	virtual bool EventPreProcessing(TVector3 &pbeam)override;
	virtual void EventPostProcessing(TVector3 &pbeam)override;
	virtual bool TrackCountTrigger(int CinC,int NinC,int CinF)override;
	virtual bool CentralFirst()override;
	virtual bool ForwardTrackProcessing(WTrack* track,TVector3 &pbeam)override;
	virtual bool CentralTrackProcessing(WTrack* track,TVector3 &pbeam)override;
private:
	vector<TH2F*> EDepHist;
	vector<TH2F*> EDepFilteredHist;
	TH1F *MissingMass;
	TH2F *MissingHist;
};
#endif