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
	auto TrackFinderFD = dynamic_cast<FDFTHTracks*>(gDataManager->GetAnalysisModule("FDFTHTracks","default"));
	if(TrackFinderFD!=0) fTrackBankFD = TrackFinderFD->GetTrackBank();
	auto CDTrackFinder = dynamic_cast<CDTracksSimple*>(gDataManager->GetAnalysisModule("CDTracksSimple","default"));
	if (CDTrackFinder!=0) fTrackBankCD = CDTrackFinder->GetTrackBank();
	He3DepKin = dynamic_cast<FDEdep2Ekin*>(gParameterManager->GetParameterObject("FDEdep2Ekin","3He"));
	fDetectorTable = dynamic_cast<CCardWDET*>(gParameterManager->GetParameterObject("CCardWDET","default")); 
}
Analysis::~Analysis(){}
void Analysis::ProcessEvent(){
	if (ProcessingCondition()){
		double event_wieght=EventWeight();
		double p_beam=PBeam();
		SpecificProcessing();
		int ChargedCountInCentral = fTrackBankCD->GetEntries(kCDC);
		int NeutralCountInCentral = fTrackBankCD->GetEntries(kCDN);
		int ChargedCountinForward = fTrackBankFD->GetEntries(kFDC);
		int NrTracks=fTrackBankFD->GetEntries();
		for (int trackindex=0; trackindex<NrTracks; trackindex++) { 
			auto track=fTrackBankFD->GetTrack(trackindex);
			auto cluster=track->GetFirstCluster();
			while (cluster) {
				//ToDo: what for?!
				if (track->Type()==kFDC && cluster->GetPlaneNumber()==kFRH2 && cluster->Element()==22) 
					return;
				if (track->Type()==kFDC && cluster->GetPlaneNumber()==kFRH3 && cluster->Element()==1) 
					return;
				cluster=track->GetNextCluster();
			}
			const int fwd_count=5; 
			double fwd_thresholds[]={0.004,0.0025,0.0025,0.0035,0.004};
			ForwardDetectorPlane fwd_planes[]={kFRH1,kFRH2,kFRH3,kFRH4,kFRH5};
			auto thresholds=[track](){
				bool res=true;
				double stop_thresholds[]={0.00018,0.00018,0.0015,0.00032,0.0003};
				ForwardDetectorPlane thr_planes[]={kFWC1,kFWC2,kFTH1,kFTH2,kFTH3};
				for(int i=0;i<5;i++)
					track->Edep(thr_planes[i])>stop_thresholds[i];
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
			if(!(thresholds()&&stop_condition))
				return;
			double EdepFRHtot=0;
			for(int i=0;i<fwd_count;i++)
				EdepFRHtot+=track->Edep(fwd_planes[i]);
			//ToDo: preselection condition for Edep
			TVector3 vec_beam;
			vec_beam.SetMagThetaPhi(p_beam,0,0);
			TLorentzVector P_b;
			P_b.SetVectM(vec_beam,m_p);
			TVector3 vec_target;
			vec_target.SetMagThetaPhi(0,0,0);
			TLorentzVector P_t;
			P_t.SetVectM(vec_target,m_d);
			TLorentzVector P_tot=P_b+P_t;
		}
	}
}

MCAnalysis::MCAnalysis():Analysis(){
	WTrackFinder *MCTrf = dynamic_cast<WTrackFinder*>(gDataManager->GetAnalysisModule("MCTrackFinder","default"));
	fMCTrackBank  = MCTrf->GetTrackBank();
	fMCVertexBank = MCTrf->GetVertexBank();
	fEventHeader = dynamic_cast<REventWmcHeader*>(gDataManager->GetDataObject("REventWmcHeader","EventHeader"));
	He3_Ekin=new TH1F("Kinetic Energy","",1000,0,2);
	He3_Theta=new TH1F("Theta","",18,0.0,30.0);
	He3_Phi=new TH1F("Phi","",36, 0,360);
	gHistoManager->Add(He3_Ekin,"E_test");
	gHistoManager->Add(He3_Theta,"Theta_test");
	gHistoManager->Add(He3_Phi,"Phi_test");
}
MCAnalysis::~MCAnalysis(){}
double MCAnalysis::EventWeight(){return fEventHeader->GetWeight();}
void MCAnalysis::SpecificProcessing(){}
bool MCAnalysis::ProcessingCondition(){
	return 
		gWasa->IsAnalysisMode(Wasa::kMCRaw)||
		gWasa->IsAnalysisMode(Wasa::kMCReco)||
		gWasa->IsAnalysisMode(Wasa::kMC);
}
double MCAnalysis::PBeam(){
	TVector3 vec_3He;
	TVector3 vec_n;
	TLorentzVector P_3He;
	TLorentzVector P_n;
	WVertexIter iterator(fMCVertexBank);
	int NrVertex=0;
	while(WVertex *vertex=dynamic_cast<WVertex*>(iterator.Next())){
		NrVertex++;
		for(int particleindex=0; particleindex<vertex->NumberOfParticles(); particleindex++){
			//ToDo: this is the calculation for another reaction
			//      one needs to implement it for my reaction
			WParticle *particle=vertex->GetParticle(particleindex);
			auto ekin=particle->GetEkin();
			auto theta=particle->GetTheta();
			auto phi=particle->GetPhi();
			if(NrVertex==1 && kHe3==particle->GetType()){
				auto p_3He=TMath::Sqrt(ekin*(ekin+2*m_3He));
				auto E_3He=TMath::Sqrt(p_3He*p_3He+m_3He*m_3He);
				vec_3He.SetMagThetaPhi(p_3He,theta,phi);
				P_3He.SetVectM(vec_3He,m_3He);
				He3_Ekin->Fill(particle->GetEkin());
				He3_Theta->Fill(particle->GetTheta()*180/3.1415926);
				He3_Phi->Fill(particle->GetPhi()*180/3.1415926);
			}
			if(NrVertex==1 && kNeutron==particle->GetType()){
				auto p_n=TMath::Sqrt(ekin*(ekin+2*m_n));
				auto E_n=TMath::Sqrt(p_n*p_n+m_n*m_n);
				vec_n.SetMagThetaPhi(p_n,theta,phi);
				P_n.SetVectM(vec_n,m_n);
			}
		}
	}
	return (vec_n+vec_3He).Mag();
}
