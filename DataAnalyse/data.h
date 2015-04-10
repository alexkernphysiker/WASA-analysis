#ifndef GXMEYPWAEFAIYNHM
#define GXMEYPWAEFAIYNHM
#include "analysis.h"
#include "math_h/interpolate.h"
class RealData:public virtual Analysis{
public:
	RealData();
	virtual ~RealData();
protected:
	virtual bool EventProcessingCondition()override;
	virtual double PBeam()override;
	virtual double EventWeight()override;
	virtual void PrepareCheck()override;
	virtual void CheckParticleTrack(ParticleType type,double Ekin,double theta, double phi)override;
private:
	REventHeader *fHeader;
	LinearInterpolation<double> p_beam;
};
#endif