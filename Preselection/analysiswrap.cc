// this file is distributed under 
// MIT license
#include <string>
#include <exception>
#include "analysiswrap.hh"
#include "config.h"
#include "log.h"
#include "he3.h"
using namespace std;
string type="";
void SetAnalysisType(string t){
	type=t;
}
Logger LOG;
AbstractAnalysis::AbstractAnalysis(){
	LOG.AddLogSubprefix("ANALYSIS APPLICATION");
	LOG.Log(LogError)<<"AnalysisWrap empty constructor. Should not be called but it is";
}
AbstractAnalysis::AbstractAnalysis(const char* name,const void*data): CAnalysisModule(name){
	Logger::SubLog log=LOG.Log(NoLog);
	if(data)
		m_data=(void*)data;
	else{
		LOG.Log(LogError)<<"Unknown analysis type";
		throw exception();
	}
}
AbstractAnalysis::~AbstractAnalysis(){
	if(m_data)
		delete (IAnalysis*)m_data;
}
void AbstractAnalysis::ProcessEvent(){
	if (fProcessed) return;
	fProcessed = kTRUE;
	((IAnalysis*)m_data)->ProcessEvent();
}
void AbstractAnalysis::Clear(Option_t *option){fProcessed=kFALSE;}
void AbstractAnalysis::Print(Option_t *option){}
void AbstractAnalysis::UserCommand(CCommand * command){}

void*GetAnalysis(){
	IAnalysis *alg=nullptr;
	if(type=="MC_He3eta")alg=new MC_He3eta();
	if((type=="MC_He3pi0")||(type=="MC_He3pi0pi0")||(type=="MC_He3pi0pi0pi0"))alg=new MC_He3pi0();
	if(type=="Data_He3")alg=new Data_He3();
	return (void*)alg;
}
void*GetReconstruction(){
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
ClassImp(ReconstructionModule);
ReconstructionModule::ReconstructionModule(const char* name):AbstractAnalysis(name,GetReconstruction()){}
ReconstructionModule::ReconstructionModule(){}
ReconstructionModule::~ReconstructionModule(){}



