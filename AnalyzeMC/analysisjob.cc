#include <Wasa.hh>
#include <CDataManager.hh>
#include <CHistoManager.hh>
#include <CParameterManager.hh>
#include <CLog.hh>
#include <CConst.hh>
#include <EmsEvent.hh>
#include <WHitBank.hh>
#include <WHitScint.hh>
#include <WVertex.hh>
#include <WVertexBank.hh>
#include <WCluster.hh>
#include <WClusterBank.hh>
#include <WClusterChamb.hh>
#include <WClusterFinder.hh>
#include <WTrack.hh>
#include <WTrackBank.hh>
#include <TString.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "analysisjob.hh"
ClassImp(AnalysisJob);
AnalysisJob::AnalysisJob(){}
AnalysisJob::AnalysisJob(const char *name):CAnalysisModule(name){
	fHeader = dynamic_cast<REventHeader*>(gDataManager->GetDataObject("REventHeader","Header"));
	TrackFinderFD = dynamic_cast<FDFTHTracks*>(gDataManager->GetAnalysisModule("FDFTHTracks","default"));
	if(TrackFinderFD) fTrackBankFD = TrackFinderFD->GetTrackBank();
	TrackFinderCD = dynamic_cast<CDTracksSimple*>(gDataManager->GetAnalysisModule("CDTracksSimple","default"));
	if (TrackFinderCD) fTrackBankCD = TrackFinderCD->GetTrackBank();
	WTrackFinder *MCTrf = dynamic_cast<WTrackFinder*>(gDataManager->GetAnalysisModule("MCTrackFinder","default"));
	fMCTrackBank  = MCTrf->GetTrackBank();fMCVertexBank = MCTrf->GetVertexBank();
	fEventHeader = dynamic_cast<REventWmcHeader*>(gDataManager->GetDataObject("REventWmcHeader","EventHeader"));
	
	He3_Ekin=new TH1F("Kinetic Energy","",1000,-0.3,0.3);
	gHistoManager->Add(He3_Ekin,"E_test");
	He3_Theta=new TH1F("Theta","",18,0.0,180.0);
	gHistoManager->Add(He3_Theta,"Theta_test");
	He3_Phi=new TH1F("Phi","",36, -180,180);
	gHistoManager->Add(He3_Phi,"Phi_test");
}
AnalysisJob::~AnalysisJob(){}
void AnalysisJob::ProcessEvent(){
	printf("marker1");
	if (fProcessed) return;
	fProcessed = kTRUE;
	if (gWasa->IsAnalysisMode(Wasa::kMCRaw)||gWasa->IsAnalysisMode(Wasa::kMCReco)||gWasa->IsAnalysisMode(Wasa::kMC)) { 
		printf("marker2\n");
		WVertexIter iterator(fMCVertexBank);
		while(WVertex *vertex=dynamic_cast<WVertex*>(iterator.Next())){
			printf("vertex:");
			for(int particleindex=0; particleindex<vertex->NumberOfParticles(); particleindex++){
				WParticle *particle=vertex->GetParticle(particleindex);
				printf("Particle: ");
				if(kHe3==particle->GetType()){
					He3_Ekin->Fill(particle->GetEkin());
					printf("3He: E = %f; Th = %f; Phi = %f",particle->GetEkin(),particle->GetTheta(), particle->GetPhi());
					He3_Theta->Fill(particle->GetTheta());
					He3_Phi->Fill(particle->GetPhi());
				}
				printf("\n");
			}
		}
	}
}
void AnalysisJob::Clear(Option_t *option){}
void AnalysisJob::Print(Option_t *option){}
void AnalysisJob::UserCommand(CCommand * command){}
