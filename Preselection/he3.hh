// this file is distributed under 
// GPL v 3.0 license
#ifndef oVonXYZj
#define oVonXYZj
#include "analysis.h"
#include "detectors.h"
#include "reconstruction.h"
class He3_production:public virtual Analysis,public ForwardDetectors{
public:
	He3_production();
	virtual ~He3_production();
protected:
	virtual bool Cuts(WTrack&&track)=0;
	virtual bool MissingMassCut(double m)=0;
	virtual bool EventPreProcessing()override;
	virtual void EventPostProcessing()override;
	virtual bool TrackCountTrigger(int CinC,int NinC,int CinF)override;
	virtual bool CentralFirst()override;
	virtual bool ForwardTrackProcessing(WTrack&&track)override;
	virtual bool CentralTrackProcessing(WTrack&&track)override;
private:
	InterpolationBasedReconstruction He3_Ekin,He3_theta,He3_phi;
	vector<TH2F*> EDepHist;
	vector<TH2F*> EDepFilteredHist;
	vector<TH1F*> MissingMassDetailed;
	TH1F *MissingMass,*DependenceOnPBeam,*P_Beam;
};

#endif 