#include "analysisjob.hh"
#include "analysis.h"
////// Calculation wraper for WASA
ClassImp(AnalysisJob);
AnalysisJob::AnalysisJob(){}
AnalysisJob::AnalysisJob(const char *name):CAnalysisModule(name){
	m_data=(void*)(new CreateAnalysis<MonteCarlo,He3eta>());
}
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
				typedef pair<WTrackBank*,function<bool(WTrack*,TVector3&)>> DetectorToProcess;
				DetectorToProcess CENTRAL=make_pair(fTrackBankCD,[](WTrack* t,TVector3& p){return CentralTrackProcessing(t,p);});
				DetectorToProcess FORWARD=make_pair(fTrackBankFD,[](WTrack* t,TVector3& p){return ForwardTrackProcessing(t,p);});
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
	TVector3 vec_3He;
	TVector3 vec_eta;
	TLorentzVector P_3He;
	TLorentzVector P_eta;
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
			if((NrVertex==2)&&(kHe3==ptype)){
				auto p_3He=TMath::Sqrt(ekin*(ekin+2*m_3He));
				auto E_3He=TMath::Sqrt(p_3He*p_3He+m_3He*m_3He);
				vec_3He.SetMagThetaPhi(p_3He,theta,phi);
				P_3He.SetVectM(vec_3He,m_3He);
			}
			if((NrVertex==2)&&(kEta==ptype)){
				auto p_eta=TMath::Sqrt(ekin*(ekin+2*m_eta));
				auto E_eta=TMath::Sqrt(p_eta*p_eta+m_eta*m_eta);
				vec_eta.SetMagThetaPhi(p_eta,theta,phi);
				P_eta.SetVectM(vec_eta,m_eta);
			}
		}
	}
	return (vec_eta+vec_3He).Mag();
}
He3eta::He3eta():Analysis(){
	FRH1vsFRH2=new TH2F("FRH1vsFRH2","",256,0,0.5,256,0,0.5);
	gHistoManager->Add(FRH1vsFRH2);
}
He3eta::~He3eta(){}
bool He3eta::EventPreProcessing(TVector3 &pbeam){
	return true;
}
void He3eta::EventPostProcessing(TVector3 &pbeam){}
bool He3eta::TrackCountTrigger(int CinC,int NinC,int CinF){
	return true;
}
bool He3eta::ForwardTrackProcessing(WTrack* track,TVector3 &pbeam){
	const int fwd_count=5; 
	double fwd_thresholds[]={0.004,0.0025,0.0025,0.0035,0.004};
	ForwardDetectorPlane fwd_planes[]={kFRH1,kFRH2,kFRH3,kFRH4,kFRH5};
	auto thresholds=[track](){
		bool res=true;
		double thresholds[]={0.00018,0.00018,0.0015,0.00032,0.0003};
		ForwardDetectorPlane thr_planes[]={kFWC1,kFWC2,kFTH1,kFTH2,kFTH3};
		for(int i=0;i<5;i++)
			track->Edep(thr_planes[i])>thresholds[i];
		return res;
	};
	auto forward_stop_plane=[track,fwd_planes,fwd_thresholds](unsigned int stopindex){
		bool res=true;
		for(int i=0;i<=stopindex;i++)
			res&=track->Edep(fwd_planes[i])>fwd_thresholds[i];
		for(int i=stopindex+1;i<fwd_count;i++)
			res&=track->Edep(fwd_planes[i])<=fwd_thresholds[i];
		return res;
	};
	bool stop_condition=false;
	for(int i=0;i<fwd_count;i++)
		stop_condition|=forward_stop_plane(i);
	if(thresholds()&&stop_condition){
		double EdepFRHtot=0;
		for(int i=0;i<fwd_count;i++)
			EdepFRHtot+=track->Edep(fwd_planes[i]);
		//ToDo: preselection condition for Edep
		
		TLorentzVector P_b;
		P_b.SetVectM(pbeam,m_p);
		TVector3 vec_target;
		vec_target.SetMagThetaPhi(0,0,0);
		TLorentzVector P_t;
		P_t.SetVectM(vec_target,m_d);
		TLorentzVector P_tot=P_b+P_t;
	}
	
	return true;
}
bool He3eta::CentralTrackProcessing(WTrack* track,TVector3 &pbeam){
	return true;
}