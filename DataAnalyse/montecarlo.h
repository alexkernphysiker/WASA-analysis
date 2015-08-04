// this file is distributed under 
// GPL v 3.0 license
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
	virtual bool GetTrueParameters(ParticleType type,double&Ekin,double&theta,double&phi)override;
private:
	WTrackBank *fMCTrackBank;
	WVertexBank *fMCVertexBank;
	REventWmcHeader *fEventHeader;
};
#endif