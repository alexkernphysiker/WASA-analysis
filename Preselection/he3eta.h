// this file is distributed under 
// GPL v 3.0 license
#ifndef EVHZDETEOYGNMFFP
#define EVHZDETEOYGNMFFP
#include "he3.hh"
class He3eta_gg_:public He3_production{
public:
	He3eta_gg_();
	virtual ~He3eta_gg_();
protected:
	virtual bool Cuts(WTrack&&track)override;
	virtual bool MissingMassCut(double m)override;
};
#endif