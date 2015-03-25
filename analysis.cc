#include "analysis.h"
using namespace std;
IAnalysis::~IAnalysis(){}
Analysis::Analysis(){
	auto TrackFinderFD = dynamic_cast<FDFTHTracks*>(gDataManager->GetAnalysisModule("FDFTHTracks","default"));
	if(TrackFinderFD!=0) fTrackBankFD = TrackFinderFD->GetTrackBank();
	auto CDTrackFinder = dynamic_cast<CDTracksSimple*>(gDataManager->GetAnalysisModule("CDTracksSimple","default"));
	if (CDTrackFinder!=0) fTrackBankCD = CDTrackFinder->GetTrackBank();
	He3DepKin = dynamic_cast<FDEdep2Ekin*>(gParameterManager->GetParameterObject("FDEdep2Ekin","3He"));
	fDetectorTable = dynamic_cast<CCardWDET*>(gParameterManager->GetParameterObject("CCardWDET","default")); 
	P_Beam=new TH1F("P_beam","",500,1.4,1.9);
	gHistoManager->Add(P_Beam,"PBeam");
}
Analysis::~Analysis(){}
typedef pair<WTrackBank*,function<bool(WTrack*,TVector3&)>> DetectorToProcess;
void Analysis::ProcessEvent(){
	if (EventProcessingCondition()){
		double event_wieght=EventWeight();
		double p_beam=PBeam();
		P_Beam->Fill(p_beam);
		TVector3 vec_beam;
		vec_beam.SetMagThetaPhi(p_beam,0,0);
		if(EventPreProcessing(vec_beam)){
			int ChargedCountInCentral = fTrackBankCD->GetEntries(kCDC);
			int NeutralCountInCentral = fTrackBankCD->GetEntries(kCDN);
			int ChargedCountinForward = fTrackBankFD->GetEntries(kFDC);
			if(TrackCountTrigger(ChargedCountInCentral,NeutralCountInCentral,ChargedCountinForward)){
				DetectorToProcess CENTRAL=make_pair(fTrackBankCD,[this](WTrack* t,TVector3& p){return CentralTrackProcessing(t,p);});
				DetectorToProcess FORWARD=make_pair(fTrackBankFD,[this](WTrack* t,TVector3& p){return ForwardTrackProcessing(t,p);});
				vector<DetectorToProcess> QUEUE;
				if(CentralFirst()){
					QUEUE.push_back(CENTRAL);
					QUEUE.push_back(FORWARD);
				}else{
					QUEUE.push_back(FORWARD);
					QUEUE.push_back(CENTRAL);
				}
				for(DetectorToProcess DETECTOR:QUEUE){
					int NrTracks=DETECTOR.first->GetEntries();
					for (int trackindex=0; trackindex<NrTracks; trackindex++) { 
						auto track=DETECTOR.first->GetTrack(trackindex);
						if(!DETECTOR.second(track,vec_beam))
							return;
					}
				}
				EventPostProcessing(vec_beam);
			}
		}
	}
}
MonteCarlo::MonteCarlo():Analysis(){
	WTrackFinder *MCTrf = dynamic_cast<WTrackFinder*>(gDataManager->GetAnalysisModule("MCTrackFinder","default"));
	fMCTrackBank  = MCTrf->GetTrackBank();
	fMCVertexBank = MCTrf->GetVertexBank();
	fEventHeader = dynamic_cast<REventWmcHeader*>(gDataManager->GetDataObject("REventWmcHeader","EventHeader"));
}
MonteCarlo::~MonteCarlo(){}
double MonteCarlo::EventWeight(){return fEventHeader->GetWeight();}
bool MonteCarlo::EventProcessingCondition(){
	return 
	gWasa->IsAnalysisMode(Wasa::kMCRaw)||
	gWasa->IsAnalysisMode(Wasa::kMCReco)||
	gWasa->IsAnalysisMode(Wasa::kMC);
}
double MonteCarlo::PBeam(){
	TVector3 result;
	WVertexIter iterator(fMCVertexBank);
	int NrVertex=0;
	while(WVertex *vertex=dynamic_cast<WVertex*>(iterator.Next())){
		NrVertex++;
		for(int particleindex=0; particleindex<vertex->NumberOfParticles(); particleindex++){
			WParticle *particle=vertex->GetParticle(particleindex);
			auto ptype=particle->GetType();
			auto ekin=particle->GetEkin();
			auto theta=particle->GetTheta();
			auto phi=particle->GetPhi();
			if(NrVertex==2)
				for(auto P:first_particles)
					if(ptype==P.first){
						double m=P.second;
						auto p=TMath::Sqrt(ekin*(ekin+2*m));
						auto E_3He=TMath::Sqrt(p*p+m*m);
						TVector3 vec;
						vec.SetMagThetaPhi(p,theta,phi);
						result=result+vec;
					}
		}
	}
	return result.Mag();
}
