#include "analysiswrap.hh"
#include "montecarlo.h"
#include "data.h"
#include "he3eta.h"
AnalysisWrap::AnalysisWrap(){}
AnalysisWrap::AnalysisWrap(const char* name): CAnalysisModule(name){}
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
/// particular types of analysis
ClassImp(MCHe3Eta);
MCHe3Eta::MCHe3Eta(){}
MCHe3Eta::MCHe3Eta(const char *name):AnalysisWrap(name){
	m_data=(void*)(new CreateAnalysis<MonteCarlo,He3eta>());
}
MCHe3Eta::~MCHe3Eta(){}
ClassImp(DataHe3Eta);
DataHe3Eta::DataHe3Eta(){}
DataHe3Eta::DataHe3Eta(const char *name):AnalysisWrap(name){
	m_data=(void*)(new CreateAnalysis<RealData,He3eta>());
}
DataHe3Eta::~DataHe3Eta(){}
