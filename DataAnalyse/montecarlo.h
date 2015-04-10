#ifndef WAFXQSAKXMIIWRNB
#define WAFXQSAKXMIIWRNB
#include "analysis.h"
class MonteCarlo:public virtual Analysis{
public:
	MonteCarlo();
	virtual ~MonteCarlo();
protected:
	virtual bool EventProcessingCondition()override;
	virtual double PBeam()override;
	virtual double EventWeight()override;
	virtual void PrepareCheck()override;
	virtual void CheckParticleTrack(ParticleType type,double Ekin,double theta, double phi)override;
private:
	struct CheckHists{
	public:
		CheckHists(ParticleType t);
		ParticleType type;
		TH1F *Ekin, *Theta, *Phi;
	};
	vector<CheckHists> check;
	WTrackBank *fMCTrackBank;
	WVertexBank *fMCVertexBank;
	REventWmcHeader *fEventHeader;
};
#endif