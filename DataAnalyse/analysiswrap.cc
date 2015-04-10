#include <string>
#include <iostream>
#include "analysiswrap.hh"
#include "montecarlo.h"
#include "data.h"
#include "he3eta.h"
using namespace std;
string type="";
void SetAnalysisType(char* t){
	type=string(t);
}
ClassImp(AnalysisWrap);
AnalysisWrap::AnalysisWrap(){}
AnalysisWrap::AnalysisWrap(const char* name): CAnalysisModule(name){
	IAnalysis *alg=nullptr;
	printf(">>>>>>>>>>>>>>>>> %s\n",type.c_str());
	if(type=="MCHe3Eta")
		alg=new CustomAnalysis<MonteCarlo,He3eta>();
	if(type=="DataHe3Eta")
		alg=new CustomAnalysis<RealData,He3eta>();
	if(alg)
		m_data=(void*)alg;
	else
		throw;
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
