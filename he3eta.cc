#include "reactions.h"
#include "detectors.h"
He3eta_gg::He3eta_gg():Analysis(){
	first_particles.push_back(make_pair(kHe3,m_3He));
	first_particles.push_back(make_pair(kEta,m_eta));
	final_particles.push_back(make_pair(kHe3,m_3He));
	final_particles.push_back(make_pair(kGamma,0));
	FRH1vsFRH2=new TH2F("FRH1vsFRH2","",256,0,0.3,256,0,0.3);
	gHistoManager->Add(FRH1vsFRH2,"FRH1vsFRH2");
	FRH2vsFRH3=new TH2F("FRH2vsFRH3","",256,0,0.3,256,0,0.3);
	gHistoManager->Add(FRH2vsFRH3,"FRH2vsFRH3");
	FRH3vsFRH4=new TH2F("FRH3vsFRH4","",256,0,0.3,256,0,0.3);
	gHistoManager->Add(FRH3vsFRH4,"FRH3vsFRH4");
	FRH4vsFRH5=new TH2F("FRH4vsFRH5","",256,0,0.3,256,0,0.3);
	gHistoManager->Add(FRH4vsFRH5,"FRH4vsFRH5");
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
	FRH1vsFRH2->Fill(track->Edep(kFRH1),track->Edep(kFRH2));
	FRH2vsFRH3->Fill(track->Edep(kFRH2),track->Edep(kFRH3));
	FRH3vsFRH4->Fill(track->Edep(kFRH3),track->Edep(kFRH4));
	FRH4vsFRH5->Fill(track->Edep(kFRH4),track->Edep(kFRH5));
	return true;
}
bool He3eta_gg::CentralTrackProcessing(WTrack* track,TVector3 &pbeam){
	return true;
}
void He3eta_gg::EventPostProcessing(TVector3 &pbeam){}