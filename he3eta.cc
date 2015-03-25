#include "reactions.h"
#include "detectors.h"
He3eta_gg::He3eta_gg():Analysis(){
	first_particles.push_back(make_pair(kHe3,m_3He));
	first_particles.push_back(make_pair(kEta,m_eta));
	final_particles.push_back(make_pair(kHe3,m_3He));
	final_particles.push_back(make_pair(kGamma,0));
	FRH1vsFRH2=new TH2F("FRH1vsFRH2","",256,0,0.5,256,0,0.5);
	gHistoManager->Add(FRH1vsFRH2,"FRH1vsFRH2");
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
	if(LastForwardDepPlane(track)>=2){
		FRH1vsFRH2->Fill(track->Edep(kFRH1),track->Edep(kFRH2));
	}
	//ToDo: preselection condition for Edep
	TLorentzVector P_b;
	P_b.SetVectM(pbeam,m_p);
	TVector3 vec_target;
	vec_target.SetMagThetaPhi(0,0,0);
	TLorentzVector P_t;
	P_t.SetVectM(vec_target,m_d);
	TLorentzVector P_tot=P_b+P_t;
	//ToDo: identify 3He and calculate missing mass
	return true;
}
bool He3eta_gg::CentralTrackProcessing(WTrack* track,TVector3 &pbeam){
	return true;
}
void He3eta_gg::EventPostProcessing(TVector3 &pbeam){}