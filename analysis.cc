#include "analysis.h"
using namespace std;
IAnalysis::~IAnalysis(){}
Analysis::Analysis(){
	checkprepared=false;
	auto TrackFinderFD = dynamic_cast<FDFTHTracks*>(gDataManager->GetAnalysisModule("FDFTHTracks","default"));
	if(TrackFinderFD!=0) fTrackBankFD = TrackFinderFD->GetTrackBank();
	auto CDTrackFinder = dynamic_cast<CDTracksSimple*>(gDataManager->GetAnalysisModule("CDTracksSimple","default"));
	if (CDTrackFinder!=0) fTrackBankCD = CDTrackFinder->GetTrackBank();
	fDetectorTable = dynamic_cast<CCardWDET*>(gParameterManager->GetParameterObject("CCardWDET","default")); 
	P_Beam=new TH1F("P_beam","",500,1.4,1.9);
	gHistoManager->Add(P_Beam,"PBeam");
}
Analysis::~Analysis(){}
void Analysis::PrepareCheck(){}
void Analysis::CheckParticleTrack(ParticleType type, double Ekin, double theta, double phi){}
typedef pair<WTrackBank*,function<bool(WTrack*,TVector3)>> DetectorToProcess;
void Analysis::ProcessEvent(){
	if(!checkprepared){
		PrepareCheck();
		checkprepared=true;
	}
	if (EventProcessingCondition()){
		double event_wieght=EventWeight();
		TVector3 p_beam;
		p_beam.SetMagThetaPhi(PBeam(),0,0);
		P_Beam->Fill(p_beam.Mag());
		if(EventPreProcessing(p_beam)){
			int ChargedCountInCentral = fTrackBankCD->GetEntries(kCDC);
			int NeutralCountInCentral = fTrackBankCD->GetEntries(kCDN);
			int ChargedCountinForward = fTrackBankFD->GetEntries(kFDC);
			if(TrackCountTrigger(ChargedCountInCentral,NeutralCountInCentral,ChargedCountinForward)){
				DetectorToProcess CENTRAL=make_pair(fTrackBankCD,[this](WTrack* t,TVector3 p){return CentralTrackProcessing(t,p);});
				DetectorToProcess FORWARD=make_pair(fTrackBankFD,[this](WTrack* t,TVector3 p){return ForwardTrackProcessing(t,p);});
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
						if(!DETECTOR.second(track,p_beam))
							return;
					}
				}
				EventPostProcessing(p_beam);
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
		if(NrVertex==2)
			for(int particleindex=0; particleindex<vertex->NumberOfParticles(); particleindex++){
				WParticle *particle=vertex->GetParticle(particleindex);
				auto ptype=particle->GetType();
				for(auto P:first_particles)
					if(ptype==P.first){
						double ekin=particle->GetEkin();
						double theta=particle->GetTheta();
						double phi=particle->GetPhi();
						double m=P.second;
						double p=sqrt(ekin*(ekin+2*m));
						TVector3 P_vec;
						P_vec.SetMagThetaPhi(p,theta,phi);
						result=result+P_vec;
					}
			}
	}
	return result.Mag();
}
MonteCarlo::CheckHists::CheckHists(ParticleType t){
	type=t;
	Ekin=new TH1F(Form("Check_E_%i",int(t)),"",500,-1,1);
	Theta=new TH1F(Form("Check_Theta_%i",int(t)),"",360,-180,180);
	Phi=new TH1F(Form("Check_Phi_%i",int(t)),"",360,-180,180);
}
void MonteCarlo::PrepareCheck(){
    Analysis::PrepareCheck();
	for(auto P:first_particles){
		CheckHists h(P.first);
		check.push_back(h);
		gHistoManager->Add(h.Ekin,h.Ekin->GetTitle());
		gHistoManager->Add(h.Theta,h.Theta->GetTitle());
		gHistoManager->Add(h.Phi,h.Phi->GetTitle());
	}
}
void MonteCarlo::CheckParticleTrack(ParticleType type, double Ekin, double theta, double phi){
	const double two_pi=2*3.1415926;
	while(phi<0)phi+=two_pi;
	while(phi>=two_pi)phi-=two_pi;
	WVertexIter iterator(fMCVertexBank);
	int NrVertex=0;
	for(CheckHists H:check)
		if(H.type==type)
			while(WVertex *vertex=dynamic_cast<WVertex*>(iterator.Next())){
				NrVertex++;
				if(NrVertex==2)
					for(int particleindex=0; particleindex<vertex->NumberOfParticles(); particleindex++){
						WParticle *particle=vertex->GetParticle(particleindex);
						auto ptype=particle->GetType();
						if(type==ptype){
							double p_ekin=particle->GetEkin();
							double p_theta=particle->GetTheta();
							double p_phi=particle->GetPhi();
							while(p_phi<0)p_phi+=two_pi;
							while(p_phi>=two_pi)p_phi-=two_pi;
							H.Ekin->Fill(Ekin-p_ekin);
							H.Theta->Fill(theta-p_theta);
							H.Phi->Fill(phi-p_phi);
						}
					}
			}
}
