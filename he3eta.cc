#include "reactions.h"
#include "detectors.h"
using namespace std;
He3eta_gg::He3eta_gg():Analysis(),ForwardDetectors(){
	first_particles.push_back(make_pair(kHe3,m_3He));
	first_particles.push_back(make_pair(kEta,m_eta));
	final_particles.push_back(make_pair(kHe3,m_3He));
	final_particles.push_back(make_pair(kGamma,0));
	int n=ForwadrPlaneCount()-1;
	for(int i=0;i<n;i++){
		string histname=ForwardPlaneName(i)+"vs"+ForwardPlaneName(i+1);
		TH2F* hist=new TH2F(histname.c_str(),"",128,0,Upper(i+1),128,0,Upper(i));
		EDepHist.push_back(hist);
		gHistoManager->Add(hist,histname.c_str());
	}
	He3DepKin = dynamic_cast<FDEdep2Ekin*>(gParameterManager->GetParameterObject("FDEdep2Ekin","3He"));
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
	int n=ForwadrPlaneCount()-1;
	for(int i=0;i<n;i++){
		EDepHist[i]->Fill(EDep(track,i+1),EDep(track,i));
	}
	int StopPlane=He3DepKin->GetStoppingPlane(track);
	if(StopPlane==kFRH1){
		int Edep2Ekin_table=0;
		double Ek3He=He3DepKin->GetEkin(kFRH1,kProton,EDep(track,kFRH1),track->Theta());
		double p3He=sqrt(Ek3He*(Ek3He+2*m_3He));
		double E3He=sqrt(p3He*p3He+m_3He*m_3He);
		double theta=track->Theta();
		double phi=track->Phi();
		CheckParticleTrack(kHe3,Ek3He,theta,phi);
		//ToDo: phi correction
		//ToDo: missing mass
	}
	return true;
}
bool He3eta_gg::CentralTrackProcessing(WTrack* track,TVector3 &pbeam){
	return true;
}
void He3eta_gg::EventPostProcessing(TVector3 &pbeam){}