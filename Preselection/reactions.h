// this file is distributed under 
// MIT license
#ifndef oVonXYZj
#define oVonXYZj
#include "analysis.h"
namespace ReactionSetup{
	enum He3ForwardModification{forData,forEta,forPi0};
	Analysis* He3_forward_analyse(He3ForwardModification);
	Analysis* He3_forward_reconstruction(He3ForwardModification);
}
#endif 