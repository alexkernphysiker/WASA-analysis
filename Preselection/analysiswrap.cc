// this file is distributed under 
// MIT license
#include <string>
#include <exception>
#include "analysiswrap.hh"
#include "config.h"
#include "log.h"
#include "analysis.h"
using namespace std;
string type="";
void SetAnalysisType(string t){type=t;}
Logger LOG;
AbstractAnalysis::AbstractAnalysis(){
	LOG.AddLogSubprefix("ANALYSIS APPLICATION");
	LOG.Log(LogError)<<"AnalysisWrap empty constructor. Should not be called but it is";
}
AbstractAnalysis::AbstractAnalysis(const char* name,const void*data):CAnalysisModule(name){
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
