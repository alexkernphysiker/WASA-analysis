// this file is distributed under 
// MIT license
#ifndef oVonXYZj
#define oVonXYZj
#include "analysis.h"
namespace ReactionSetup{
	enum He3Modification{forData,forEta,forPi0};
	Analysis* He3_forward_analyse(He3Modification);
	Analysis* He3_forward_reconstruction(He3Modification);
}
#endif 