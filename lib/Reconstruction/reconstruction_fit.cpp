// this file is distributed under 
// MIT license
#include <fstream>
#include <math_h/exception_math_h.h>
#include <str_get.h>
#include "reconstruction_fit.h"
namespace SimulationDataProcess{
	using namespace std;
	using namespace Genetic;
	string SimulationDataPath(){
		static string str=ENV(PRESEL_DATA)+string("/../Reconstruction/");
		return str;
	}
};
