// this file is distributed under 
// MIT license
#include "analysismodule.hh"
#include "he3.h"
extern string type;
void*GetAnalysis(){
	IAnalysis *alg=nullptr;
	if(type=="MC_He3eta")alg=new MC_He3eta();
	if((type=="MC_He3pi0")||(type=="MC_He3pi0pi0")||(type=="MC_He3pi0pi0pi0"))alg=new MC_He3pi0();
	if(type=="Data_He3")alg=new Data_He3();
	return (void*)alg;
}
ClassImp(AnalysisModule);
AnalysisModule::AnalysisModule(){}
AnalysisModule::AnalysisModule(const char* name): AbstractAnalysis(name,GetAnalysis()){}
AnalysisModule::~AnalysisModule(){}