#include <Wasa.hh>
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
#include <CAnalysisModule.hh>
#include <FDEdep2Ekin.hh>
#include <CCardWDET.hh>
#include "analysis.h"
const Double_t m_p=0.938272;//[GeV]
const Double_t m_n=0.93956;//[GeV]
const Double_t m_d=1.875613;//[GeV]
const Double_t m_3He=2.808950;//[GeV]
const Double_t m_4He=3.7264225;//[GeV]
const Double_t m_eta=0.547853;//[GeV]
enum TrackType{kFDN=1,kFDC=2,kCDN=11,kCDC=12};//FD neutral, FD charged, CD neutral, CD charged
enum ParticleType{kDummy=0,kGamma=1,kElectron=2,kPositron=3,kPi0=7,kPiPlus=8,kPiMinus=9,kNeutron=13,kProton=14,kDeuteron=45,kTriton=46,kHe3=49};
enum ForwardDetectorPlane{kFWC1=10,kFWC2=11,kFTH1=1,kFTH2=2,kFTH3=3,kFRH1=4,kFRH2=5,kFRH3=6,kFRH4=7,kFRH5=8,kFVH=9};
Analysis::Analysis(){
	He3_Ekin=new TH1F("Kinetic Energy","",1000,0,2);
	He3_Theta=new TH1F("Theta","",18,0.0,30.0);
	He3_Phi=new TH1F("Phi","",36, 0,360);
}
Analysis::~Analysis(){}
void Analysis::Init(CDataManager *dataManager, CHistoManager *histoManager){
	fHeader = dynamic_cast<REventHeader*>(dataManager->GetDataObject("REventHeader","Header"));
	TrackFinderFD = dynamic_cast<FDFTHTracks*>(dataManager->GetAnalysisModule("FDFTHTracks","default"));
	if(TrackFinderFD) fTrackBankFD = TrackFinderFD->GetTrackBank();
	TrackFinderCD = dynamic_cast<CDTracksSimple*>(dataManager->GetAnalysisModule("CDTracksSimple","default"));
	if (TrackFinderCD) fTrackBankCD = TrackFinderCD->GetTrackBank();
	WTrackFinder *MCTrf = dynamic_cast<WTrackFinder*>(dataManager->GetAnalysisModule("MCTrackFinder","default"));
	fMCTrackBank  = MCTrf->GetTrackBank();fMCVertexBank = MCTrf->GetVertexBank();
	fEventHeader = dynamic_cast<REventWmcHeader*>(dataManager->GetDataObject("REventWmcHeader","EventHeader"));
	histoManager->Add(He3_Ekin,"E_test");
	histoManager->Add(He3_Theta,"Theta_test");
	histoManager->Add(He3_Phi,"Phi_test");
}
void Analysis::Processevent(Wasa *wasa){
	if (wasa->IsAnalysisMode(Wasa::kMCRaw)||wasa->IsAnalysisMode(Wasa::kMCReco)||wasa->IsAnalysisMode(Wasa::kMC)) {
		/*WVertexIter iterator(vertexBank);
		while(WVertex *vertex=dynamic_cast<WVertex*>(iterator.Next())){
			for(int particleindex=0; particleindex<vertex->NumberOfParticles(); particleindex++){
				WParticle *particle=vertex->GetParticle(particleindex);
				if(kHe3==particle->GetType()){
					He3_Ekin->Fill(particle->GetEkin());
					He3_Theta->Fill(particle->GetTheta()*180/3.1415926);
					He3_Phi->Fill(particle->GetPhi()*180/3.1415926);
				}
			}
		}*/
	}
}
