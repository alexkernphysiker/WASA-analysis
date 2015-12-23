// this file is distributed under 
// MIT license
#include "reconstructionmodule.hh"
#include "he3.h"
#include "config.h"
string type="";
void SetReconstructionType(string t){type=t;}
void*GetReconstruction(){
	IAnalysis *alg=nullptr;
	if(type=="MC_He3eta")alg=new MC_He3eta();
	if((type=="MC_He3pi0")||(type=="MC_He3pi0pi0")||(type=="MC_He3pi0pi0pi0"))alg=new MC_He3pi0();
	if(type=="Data_He3")alg=new Data_He3();
	return (void*)alg;
}
ClassImp(ReconstructionModule);
ReconstructionModule::ReconstructionModule(const char* name):AbstractAnalysis(name,GetReconstruction()){}
ReconstructionModule::ReconstructionModule(){}
ReconstructionModule::~ReconstructionModule(){}