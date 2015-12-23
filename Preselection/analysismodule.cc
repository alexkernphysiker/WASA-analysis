// this file is distributed under 
// MIT license
#include "analysismodule.hh"
#include "config.h"
using namespace std;
Logger LOG;
IAnalysis*gAnalysis=nullptr;
void SetAnalysis(IAnalysis*analysis){gAnalysis=analysis;}
ClassImp(AnalysisModule);
AnalysisModule::AnalysisModule(){
	LOG.AddLogSubprefix("ANALYSIS APPLICATION");
	LOG.Log(LogError)<<"AnalysisWrap empty constructor. Should not be called but it is";
}
AnalysisModule::AnalysisModule(const char* name):CAnalysisModule(name){
	Logger::SubLog log=LOG.Log(NoLog);
	if(gAnalysis)
		m_data=(void*)gAnalysis;
	else{
		LOG.Log(LogError)<<"Unknown analysis type";
		throw exception();
	}
}
AnalysisModule::~AnalysisModule(){
	if(m_data)
		delete (IAnalysis*)m_data;
}
void AnalysisModule::ProcessEvent(){
	if (fProcessed) return;
	fProcessed = kTRUE;
	((IAnalysis*)m_data)->ProcessEvent();
}
void AnalysisModule::Clear(Option_t *option){fProcessed=kFALSE;}
void AnalysisModule::Print(Option_t *option){}
void AnalysisModule::UserCommand(CCommand * command){}
