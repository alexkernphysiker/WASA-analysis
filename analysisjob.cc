#include <memory>
#include "analysisjob.hh"
#include <TLorentzVector.h>
#include <CAnalysisModule.hh>
#include <CDataManager.hh>
#include <FDEdep2Ekin.hh>
#include <CCardWDET.hh>
#include "math_h/sigma.h"
using namespace std;
const Double_t m_p=0.938272;//[GeV]
const Double_t m_n=0.93956;//[GeV]
const Double_t m_d=1.875613;//[GeV]
const Double_t m_3He=2.808950;//[GeV]
const Double_t m_4He=3.7264225;//[GeV]
const Double_t m_eta=0.547853;//[GeV]
enum TrackType{kFDN=1,kFDC=2,kCDN=11,kCDC=12};//FD neutral, FD charged, CD neutral, CD charged
enum ParticleType{kDummy=0,kGamma=1,kElectron=2,kPositron=3,kPi0=7,kPiPlus=8,kPiMinus=9,kNeutron=13,kProton=14,kDeuteron=45,kTriton=46,kHe3=49};
enum ForwardDetectorPlane{kFWC1=10,kFWC2=11,kFTH1=1,kFTH2=2,kFTH3=3,kFRH1=4,kFRH2=5,kFRH3=6,kFRH4=7,kFRH5=8,kFVH=9};
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
	He3_Ekin=new TH1F("Kinetic Energy","",1000,0,2);
	He3_Theta=new TH1F("Theta","",18,0.0,30.0);
	He3_Phi=new TH1F("Phi","",36, 0,360);
	gHistoManager->Add(He3_Ekin,"E_test");
	gHistoManager->Add(He3_Theta,"Theta_test");
	gHistoManager->Add(He3_Phi,"Phi_test");

	auto sig=make_shared<Sigma<double>>();
}
AnalysisJob::~AnalysisJob(){
}
void AnalysisJob::ProcessEvent(){
	if (fProcessed) return;
	fProcessed = kTRUE;
	/*if (gWasa->IsAnalysisMode(Wasa::kMCRaw)||gWasa->IsAnalysisMode(Wasa::kMCReco)||gWasa->IsAnalysisMode(Wasa::kMC)) {
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
	}*/
}
void AnalysisJob::Clear(Option_t *option){
	fProcessed=kFALSE;
}
void AnalysisJob::Print(Option_t *option){}
void AnalysisJob::UserCommand(CCommand * command){}
