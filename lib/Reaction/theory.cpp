// this file is distributed under 
// MIT license
#include <unistd.h>
#include "../../phys_constants.h"
#include "theory.h"
namespace Reactions{
	double sigmaHe3eta(double p_beam){
		if(p_beam<=p_he3_eta_threshold)return 0;
		else return 400;
	}
};