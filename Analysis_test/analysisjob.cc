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
#include<TVector3.h>
#include<TLorentzVector.h>
#include "analysisjob.hh"
ClassImp(AnalysisJob);
AnalysisJob::AnalysisJob(){}
AnalysisJob::AnalysisJob(const char *name):CAnalysisModule(name){
	fHeader = dynamic_cast<REventHeader*>(gDataManager->GetDataObject("REventHeader","Header"));
	TrackFinderFD = dynamic_cast<FDFTHTracks*>(gDataManager->GetAnalysisModule("FDFTHTracks","default"));
	if(TrackFinderFD) fTrackBankFD = TrackFinderFD->GetTrackBank();
	TrackFinderCD = dynamic_cast<CDTracksSimple*>(gDataManager->GetAnalysisModule("CDTracksSimple","default"));
	if (TrackFinderCD) fTrackBankCD = TrackFinderCD->GetTrackBank();
	fFDEdep2Ekin = dynamic_cast<FDEdep2Ekin*>(gParameterManager->GetParameterObject("FDEdep2Ekin","3He"));
	WTrackFinder *MCTrf = dynamic_cast<WTrackFinder*>(gDataManager->GetAnalysisModule("MCTrackFinder","default"));
	fMCTrackBank  = MCTrf->GetTrackBank();fMCVertexBank = MCTrf->GetVertexBank();
	fEventHeader = dynamic_cast<REventWmcHeader*>(gDataManager->GetDataObject("REventWmcHeader","EventHeader"));
	//TODO setup histograms
}
AnalysisJob::~AnalysisJob(){}
void AnalysisJob::ProcessEvent(){
	if (fProcessed) return;
	fProcessed = kTRUE;
	Double_t weight=1;
	if (gWasa->IsAnalysisMode(Wasa::kMCRaw)||gWasa->IsAnalysisMode(Wasa::kMCReco)||gWasa->IsAnalysisMode(Wasa::kMC))
		weight=fEventHeader->GetWeight();
	//TODO get event parameters and analyse them
}
void AnalysisJob::Clear(Option_t *option){}
void AnalysisJob::Print(Option_t *option){}
void AnalysisJob::UserCommand(CCommand * command){}
