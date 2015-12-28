// this file is distributed under 
// MIT license
#include <fstream>
#include <str_get.h>
#include "read_simulation.h"
using namespace std;
using namespace Genetic;
string SimulationDataPath(){
	static string str="";
	if(str=="")str=ENV(PRESEL_DATA)+string("/../Reconstruction/");
	return str;
}
