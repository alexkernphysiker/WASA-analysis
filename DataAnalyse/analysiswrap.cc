#include "analysiswrap.hh"
#include "montecarlo.h"
#include "data.h"
#include "he3eta.h"
ClassImp(MCHe3Eta);
MCHe3Eta::MCHe3Eta(){}
MCHe3Eta::MCHe3Eta(const char *name):CAnalysisModule(name){
	m_data=(void*)(new CreateAnalysis<MonteCarlo,He3eta>());
}
MCHe3Eta::~MCHe3Eta(){delete (IAnalysis*)m_data;}
void MCHe3Eta::ProcessEvent(){
	if (fProcessed) return;
	fProcessed = kTRUE;
	((IAnalysis*)m_data)->ProcessEvent();
}
void MCHe3Eta::Clear(Option_t *option){fProcessed=kFALSE;}
void MCHe3Eta::Print(Option_t *option){}
void MCHe3Eta::UserCommand(CCommand * command){}
