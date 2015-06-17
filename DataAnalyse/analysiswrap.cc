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
	Logger::SubLog log=LOG.getSubLog("AnalysisWrap empty constructor");
	log.Message(LogWarning,"Should not appear but appears");
}
AnalysisWrap::AnalysisWrap(const char* name): CAnalysisModule(name){
	Logger::SubLog log=LOG.getSubLog("AnalysisWrap normal constructor");
	IAnalysis *alg=nullptr;
	log.Message(NoLog,"Analysis type:");
	log.Message(NoLog,type);
	if(type=="MC_He3eta_gg_")
		alg=new CustomAnalysis<MonteCarlo,He3eta_gg_>();
	if(type=="Data_He3eta_gg_")
		alg=new CustomAnalysis<RealData,He3eta_gg_>();
	if(alg)
		m_data=(void*)alg;
	else{
		log.Message(LogError,"Unknown type");
		throw exception();
	}
}
AnalysisWrap::~AnalysisWrap(){
	Logger::SubLog log=LOG.getSubLog("AnalysisWrap destructor");
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
