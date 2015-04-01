#include "reactions.h"
#include "detectors.h"
using namespace std;
He3eta_gg::He3eta_gg():Analysis(),ForwardDetectorRoutines<kProton,kFRH1>("3He"){
	first_particles.push_back(make_pair(kHe3,m_3He));
	first_particles.push_back(make_pair(kEta,m_eta));
	final_particles.push_back(make_pair(kHe3,m_3He));
	final_particles.push_back(make_pair(kGamma,0));
	for(int i=0;i<(ForwadrPlaneCount()-1);i++){
		string histname=ForwardPlaneName(i)+"vs"+ForwardPlaneName(i+1);
		TH2F* hist=new TH2F(histname.c_str(),"",128,0,Upper(i+1),128,0,Upper(i));
		EDepHist.push_back(hist);
		gHistoManager->Add(hist,histname.c_str());
		histname=histname+"_filtered";
		hist=new TH2F(histname.c_str(),"",128,0,Upper(i+1),128,0,Upper(i));
		EDepFilteredHist.push_back(hist);
		gHistoManager->Add(hist,histname.c_str());
	}
	SetGettableFunction([](int& StopPlane,int& Edep2Ekin_table){
		//Jakieś dziwne są te warunki
		int last_plane=StopPlane;
		if(StopPlane==9 || StopPlane==8) {
			Edep2Ekin_table=17;
			last_plane=7;
		}else{
			Edep2Ekin_table=StopPlane;
			last_plane=StopPlane;
		}
		if(Edep2Ekin_table==3) 
			return false;
		else{
			StopPlane=last_plane;
			return true;
		}
	});
	SetCorrectionCoefficient(1.003);//po co?
}
He3eta_gg::~He3eta_gg(){}
bool He3eta_gg::EventPreProcessing(TVector3 &pbeam){
	return true;
}
bool He3eta_gg::TrackCountTrigger(int CinC,int NinC,int CinF){
	return true;
}
bool He3eta_gg::CentralFirst(){
	return false;
}
bool He3eta_gg::ForwardTrackProcessing(WTrack* track,TVector3 &pbeam){
	for(int i=0;i<(ForwadrPlaneCount()-1);i++)
		EDepHist[i]->Fill(EDep(track,i+1),EDep(track,i));
	auto StopPlane=StoppingPlane(track);
	if((StopPlane==kFRH1)&&(track->Edep(kFTH3)>0.01)){
		for(int i=0;i<(ForwadrPlaneCount()-1);i++)
			EDepFilteredHist[i]->Fill(EDep(track,i+1),EDep(track,i));
		double Ek3He=0;
		if(ReconstructEkin(track,Ek3He)){
			double theta=track->Theta();
			double phi=track->Phi();
			CheckParticleTrack(kHe3,Ek3He,theta,phi);
			//ToDo: phi correction
			double p3He=sqrt(Ek3He*(Ek3He+2*m_3He));
			double E3He=sqrt(p3He*p3He+m_3He*m_3He);
			//ToDo: missing mass
		}
	}
	return true;
}
bool He3eta_gg::CentralTrackProcessing(WTrack* track,TVector3 &pbeam){
	return true;
}
void He3eta_gg::EventPostProcessing(TVector3 &pbeam){}