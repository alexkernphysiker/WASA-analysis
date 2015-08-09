// this file is distributed under 
// GPL v 3.0 license
#include "he3eta.h"
He3eta_gg_::He3eta_gg_(){}
He3eta_gg_::~He3eta_gg_(){}
bool He3eta_gg_::Cuts(WTrack&& track){
	return true
	&&(EDep(static_right(track),kFRH1)>0.08)
	&&(EDep(static_right(track),kFRH1)<0.23)
	;
}
bool He3eta_gg_::MissingMassCut(double m){
	return true;
}
