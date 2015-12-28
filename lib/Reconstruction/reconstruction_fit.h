// this file is distributed under 
// MIT license
#ifndef JUDIQVAJ
# define JUDIQVAJ
#include <iostream>
#include <string>
#include <sstream>
#include <Genetic/fit.h>
#include "../../reconstruction_types.h"
template<class FITFUNC>
void ProcessFit(std::string reconstructionname, std::shared_ptr<Genetic::IInitialConditions>init,std::shared_ptr<Genetic::IParamCheck>filter){
	using namespace std;
	using namespace Genetic;
	auto points=make_shared<FitPoints>();
	ifstream file;
	file.open(reconstructionname+".simulation.txt");
	if(file){
		string line;
		while(getline(file,line)){
			istringstream str(line);
			ParamSet X;
			str>>X;
			double y;
			X>>y;
			points<<Point(X,y);
		}
		file.close();
	}
	
}
#endif