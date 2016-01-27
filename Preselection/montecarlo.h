// this file is distributed under 
// MIT license
#ifndef WAFXQSAKXMIIWRNB
#define WAFXQSAKXMIIWRNB
#include "analysis.h"
class MonteCarlo:public Analysis{
public:
	MonteCarlo();
	virtual ~MonteCarlo();
protected:
	virtual bool DataTypeSpecificEventAnalysis()override;
	virtual bool DataSpecificTriggerCheck(int n)const override;
private:
	REventHeader *fHeader;
	REventWmcHeader *fEventHeader;
	WTrackBank *fMCTrackBank;
	WVertexBank *fMCVertexBank;
};
#endif