// this file is distributed under 
// GPL v 3.0 license
#include <unistd.h>
#include "../../phys_constants.h"
#include "theory.h"
double sigmaHe3eta(double p_beam){
	if(p_beam<=p_he3_eta_threshold)return 0;
	else return 400;
}