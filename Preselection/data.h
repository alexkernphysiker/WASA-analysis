// this file is distributed under 
// GPL v 3.0 license
#ifndef GXMEYPWAEFAIYNHM
#define GXMEYPWAEFAIYNHM
#include "analysis.h"
#include "math_h/interpolate.h"
class RealData:public virtual Analysis{
public:
	RealData();
	virtual ~RealData();
protected:
    virtual void PrepairForEventAnalysis()override;
	virtual bool EventProcessingCondition()override;
private:
	REventHeader *fHeader;
	LinearInterpolation<double> p_beam;
};
#endif