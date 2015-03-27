#ifndef EVHZDETEOYGNMFFP
#define EVHZDETEOYGNMFFP
#include "analysis.h"
#include "detectors.h"
const double m_p=0.938272;//[GeV]
const double m_n=0.93956;//[GeV]
const double m_d=1.875613;//[GeV]
const double m_3He=2.808950;//[GeV]
const double m_eta=0.547853;//[GeV]
class He3eta_gg:public virtual ForwardDetectors{
public:
	He3eta_gg();
	virtual ~He3eta_gg();
protected:
	virtual bool EventPreProcessing(TVector3 &pbeam)override;
	virtual void EventPostProcessing(TVector3 &pbeam)override;
	virtual bool TrackCountTrigger(int CinC,int NinC,int CinF)override;
	virtual bool CentralFirst()override;
	virtual bool ForwardTrackProcessing(WTrack* track,TVector3 &pbeam)override;
	virtual bool CentralTrackProcessing(WTrack* track,TVector3 &pbeam)override;
private:
	vector<TH2F*> EDepHist;
};
#endif