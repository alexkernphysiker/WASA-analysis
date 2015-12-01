// this file is distributed under 
// MIT license
#ifndef WAFXQSAKXMIIWRNB
#define WAFXQSAKXMIIWRNB
#include "analysis.h"
class MonteCarlo:public virtual Analysis{
public:
	MonteCarlo();
	virtual ~MonteCarlo();
protected:
	virtual void PrepairForEventAnalysis()override;
	virtual bool EventProcessingCondition()override;
private:
	WTrackBank *fMCTrackBank;
	WVertexBank *fMCVertexBank;
	REventWmcHeader *fEventHeader;
};
#endif