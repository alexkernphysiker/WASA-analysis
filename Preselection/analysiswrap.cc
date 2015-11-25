// this file is distributed under 
// GPL v 3.0 license
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
ClassImp(AnalysisWrap);
AnalysisWrap::AnalysisWrap(){
	LOG.AddLogSubprefix("ANALYSIS APPLICATION");
	LOG.Log(LogError)<<"AnalysisWrap empty constructor. Should not be called but it is";
}
AnalysisWrap::AnalysisWrap(const char* name): CAnalysisModule(name){
	Logger::SubLog log=LOG.Log(NoLog);
	IAnalysis *alg=nullptr;
	log<<"Analysis type:"<<type;
	if(type=="MC_He3eta")alg=new MC_He3eta();
	if((type=="MC_He3pi0")||(type=="MC_He3pi0pi0")||(type=="MC_He3pi0pi0pi0"))alg=new MC_He3pi0();
	if(type=="Data_He3")alg=new Data_He3();
	if(alg)
		m_data=(void*)alg;
	else{
		LOG.Log(LogError)<<"Unknown analysis type";
		throw exception();
	}
}
AnalysisWrap::~AnalysisWrap(){
	if(m_data)
		delete (IAnalysis*)m_data;
}
void AnalysisWrap::ProcessEvent(){
	if (fProcessed) return;
	fProcessed = kTRUE;
	((IAnalysis*)m_data)->ProcessEvent();
}
void AnalysisWrap::Clear(Option_t *option){fProcessed=kFALSE;}
void AnalysisWrap::Print(Option_t *option){}
void AnalysisWrap::UserCommand(CCommand * command){}
