#include "analysisjob.hh"
#include "analysis.h"
#include "reactions.h"
ClassImp(AnalysisJob);
AnalysisJob::AnalysisJob(){}
AnalysisJob::AnalysisJob(const char *name):CAnalysisModule(name){
	m_data=(void*)(new CreateAnalysis<MonteCarlo,He3eta_gg>());
}
AnalysisJob::~AnalysisJob(){delete (IAnalysis*)m_data;}
void AnalysisJob::ProcessEvent(){
	if (fProcessed) return;
	fProcessed = kTRUE;
	((IAnalysis*)m_data)->ProcessEvent();
}
void AnalysisJob::Clear(Option_t *option){fProcessed=kFALSE;}
void AnalysisJob::Print(Option_t *option){}
void AnalysisJob::UserCommand(CCommand * command){}
