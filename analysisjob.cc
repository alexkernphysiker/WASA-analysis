#include "analysisjob.hh"
#include "AnalyseBase/wrap.h"
ClassImp(AnalysisJob);
AnalysisJob::AnalysisJob(){}
AnalysisJob::AnalysisJob(const char *name):CAnalysisModule(name){
	job=CreateAnalysis();
	job->Init(gDataManager,gHistoManager);
}
AnalysisJob::~AnalysisJob(){delete job;}
void AnalysisJob::ProcessEvent(){
	if (fProcessed) return;
	fProcessed = kTRUE;
	job->Processevent(gWasa);
}
void AnalysisJob::Clear(Option_t *option){
	fProcessed=kFALSE;
}
void AnalysisJob::Print(Option_t *option){}
void AnalysisJob::UserCommand(CCommand * command){}
