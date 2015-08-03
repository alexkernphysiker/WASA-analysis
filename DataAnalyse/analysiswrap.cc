// this file is distributed under 
// GPL v 3.0 license
#include <string>
#include <exception>
#include "config.h"
#include "log.h"
#include "analysiswrap.hh"
#include "montecarlo.h"
#include "data.h"
#include "he3eta.h"
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
	if(type=="MC_He3eta_gg_")
		alg=new CustomAnalysis<MonteCarlo,He3eta_gg_>();
	if(type=="Data_He3eta_gg_")
		alg=new CustomAnalysis<RealData,He3eta_gg_>();
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
