// this file is distributed under 
// MIT license
#include <functional>
#include "analysismodule.hh"
#include "config.h"
#include "he3.h"
using namespace std;
Logger LOG;
string type;
void SetAnalysisType(std::string t){
	type=t;
}

ClassImp(AnalysisModule);
AnalysisModule::AnalysisModule(){
	LOG.AddLogSubprefix("ANALYSIS APPLICATION");
	LOG.Log(LogError)<<"AnalysisWrap empty constructor. Should not be called but it is";
}
AnalysisModule::AnalysisModule(const char* name):CAnalysisModule(name){
	Logger::SubLog log=LOG.Log(NoLog);
	if("RE_He3eta"==type)
		m_data= new RE_He3eta();
	if(
		("RE_He3pi0"==type)||
		("RE_He3pi0pi0"==type)||
		("RE_He3pi0pi0pi0"==type)
	)
		m_data=new RE_He3pi0();
	if("MC_He3eta"==type)
		m_data=new MC_He3eta();
	if(
		("MC_He3pi0"==type)||
		("MC_He3pi0pi0"==type)||
		("MC_He3pi0pi0pi0"==type)
	)
		m_data=new MC_He3pi0();
	if("Data_He3"==type)
		m_data=new Data_He3();
}
AnalysisModule::~AnalysisModule(){
	if(m_data)
		delete m_data;
}
void AnalysisModule::ProcessEvent(){
	if (fProcessed) return;
	fProcessed = kTRUE;
	((IAnalysis*)m_data)->ProcessEvent();
}
void AnalysisModule::Clear(Option_t *option){fProcessed=kFALSE;}
void AnalysisModule::Print(Option_t *option){}
void AnalysisModule::UserCommand(CCommand * command){}
