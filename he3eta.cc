#include <iostream>
#include "reactions.h"
#include "detectors.h"
using namespace std;
He3eta::He3eta():Analysis(),ForwardDetectorRoutines("3He"){
	first_particles.push_back(make_pair(kHe3,m_3He));
	first_particles.push_back(make_pair(kEta,m_eta));
	final_particles.push_back(make_pair(kHe3,m_3He));
	final_particles.push_back(make_pair(kGamma,0));
	for(int i=0;i<(ForwadrPlaneCount()-1);i++){
		string histname=ForwardPlaneName(i)+"vs"+ForwardPlaneName(i+1);
		TH2F* hist;
		auto makehist=[this,i,&hist,&histname](vector<TH2F*> &vect){
			hist=new TH2F(histname.c_str(),"",128,0,UpperByIndex(i+1),128,0,UpperByIndex(i));
			vect.push_back(hist);
			gHistoManager->Add(hist,histname.c_str());
		};
		makehist(EDepHist);
		histname=histname+"_filtered";
		makehist(EDepFilteredHist);
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
	MissingMass=new TH1F("MissingMass","",500,0.4,0.6);
	gHistoManager->Add(MissingMass,"MissingMass");
	MissingHist=new TH2F("MissingMass_vs_Energy","",250,0.4,0.6,500,0.0,1.0);
	gHistoManager->Add(MissingHist,"MissingMass_vs_Energy");
}
He3eta::~He3eta(){}
bool He3eta::EventPreProcessing(TVector3 &pbeam){
	return true;
}
bool He3eta::TrackCountTrigger(int CinC,int NinC,int CinF){
	return true;
}
bool He3eta::CentralFirst(){
	return false;
}
bool He3eta::ForwardTrackProcessing(WTrack* track,TVector3 &p_beam){
	for(int i=0;i<(ForwadrPlaneCount()-1);i++)
		EDepHist[i]->Fill(EDep(track,i+1),EDep(track,i));
	if(
		(StoppingPlane(track)==kFRH1)
		&&UpperThresholdUpTo(track,kFRH1)
		&&(track->Edep(kFTH3)>0.01)
	){
		for(int i=0;i<(ForwadrPlaneCount()-1);i++)
			EDepFilteredHist[i]->Fill(EDep(track,i+1),EDep(track,i));
		double Ek=0;
		if(ReconstructEkin(track,Ek)){
			double theta=track->Theta();
			double phi=track->Phi();
			//ToDo: phi correction
			CheckParticleTrack(kHe3,Ek,theta,phi);
			double p3He=sqrt(Ek*(Ek+2*m_3He));
			double E3He=sqrt(p3He*p3He+m_3He*m_3He);
			TVector3 p_He3;
			p_He3.SetMagThetaPhi(p3He,theta,phi);
			TLorentzVector P_He3;
			P_He3.SetVectM(p_He3,m_3He);
			TLorentzVector P_Total;
			{
				TLorentzVector P_Beam;
				TLorentzVector P_Target;
				P_Beam.SetVectM(p_beam,m_p);
				TVector3 ptarget;
				ptarget.SetMagThetaPhi(0,0,0);
				P_Target.SetVectM(ptarget,m_d);
				P_Total=P_Beam+P_Target;
			}
			TLorentzVector P_Missing=P_Total-P_He3;
			MissingMass->Fill(P_Missing.M());
			MissingHist->Fill(P_Missing.M(),P_Missing.E());
		}
	}
	return true;
}
bool He3eta::CentralTrackProcessing(WTrack* track,TVector3 &pbeam){
	return true;
}
void He3eta::EventPostProcessing(TVector3 &pbeam){}
