#include "analysisjob.hh"
#include "analysis.h"
////// Calculation wraper for WASA
ClassImp(AnalysisJob);
AnalysisJob::AnalysisJob(){}
AnalysisJob::AnalysisJob(const char *name):CAnalysisModule(name){m_data=(void*)(new Analysis());}
AnalysisJob::~AnalysisJob(){delete (IAnalysis*)m_data;}
void AnalysisJob::ProcessEvent(){
	if (fProcessed) return;
	fProcessed = kTRUE;
	((IAnalysis*)m_data)->ProcessEvent();
}
void AnalysisJob::Clear(Option_t *option){fProcessed=kFALSE;}
void AnalysisJob::Print(Option_t *option){}
void AnalysisJob::UserCommand(CCommand * command){}
IAnalysis::~IAnalysis(){}
using namespace std;
///// Calculation implementation
Analysis::Analysis(){
	WTrackFinder *MCTrf = dynamic_cast<WTrackFinder*>(gDataManager->GetAnalysisModule("MCTrackFinder","default"));
	fMCTrackBank  = MCTrf->GetTrackBank();
	fMCVertexBank = MCTrf->GetVertexBank();
	He3_Ekin=new TH1F("Kinetic Energy","",1000,0,2);
	He3_Theta=new TH1F("Theta","",18,0.0,30.0);
	He3_Phi=new TH1F("Phi","",36, 0,360);
	gHistoManager->Add(He3_Ekin,"E_test");
	gHistoManager->Add(He3_Theta,"Theta_test");
	gHistoManager->Add(He3_Phi,"Phi_test");
}
Analysis::~Analysis(){}
void Analysis::ProcessEvent(){
	if (gWasa->IsAnalysisMode(Wasa::kMCRaw)||gWasa->IsAnalysisMode(Wasa::kMCReco)||gWasa->IsAnalysisMode(Wasa::kMC)) {
		WVertexIter iterator(fMCVertexBank);
		while(WVertex *vertex=dynamic_cast<WVertex*>(iterator.Next())){
			for(int particleindex=0; particleindex<vertex->NumberOfParticles(); particleindex++){
				WParticle *particle=vertex->GetParticle(particleindex);
				if(kHe3==particle->GetType()){
					He3_Ekin->Fill(particle->GetEkin());
					He3_Theta->Fill(particle->GetTheta()*180/3.1415926);
					He3_Phi->Fill(particle->GetPhi()*180/3.1415926);
				}
			}
		}
	}
}
