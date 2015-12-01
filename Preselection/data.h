// this file is distributed under 
// MIT license
#ifndef GXMEYPWAEFAIYNHM
#define GXMEYPWAEFAIYNHM
#include "analysis.h"
#include "reconstruction.h"
class RealData:public virtual Analysis{
public:
	RealData();
	virtual ~RealData();
protected:
    virtual void PrepairForEventAnalysis()override;
	virtual bool EventProcessingCondition()override;
private:
	REventHeader *fHeader;
	InterpolationBasedReconstruction<> BeamMomenta;
};
#endif