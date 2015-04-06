#include <iostream>
#include "he3eta.h"
#include "detectors.h"
#include "../General/phys_constants.h"
using namespace std;
He3eta::He3eta():Analysis(),ForwardDetectorRoutines("3He"){
	first_particles.push_back(make_pair(kHe3,m_3He));
	first_particles.push_back(make_pair(kEta,m_eta));
	final_particles.push_back(make_pair(kHe3,m_3He));
	final_particles.push_back(make_pair(kGamma,0));
	for(int i=0,n=ForwadrPlaneCount()-1;i<n;i++){
		string histname=ForwardPlaneName(i)+"_vs_"+ForwardPlaneName(i+1);
		TH2F* hist;
		auto makehist=[this,i,&hist,&histname](vector<TH2F*> &vect){
			hist=new TH2F(histname.c_str(),"",128,0,UpperByIndex(i+1),128,0,UpperByIndex(i));
			vect.push_back(hist);
			gHistoManager->Add(hist,"E_dep");
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
	DependenceOnPBeam=new TH1F("DependenceOnBeam","",P_hist_param);
	gHistoManager->Add(DependenceOnPBeam,"OutputEventsCount");
#define missingmassparam 600,0.3,0.6
	MissingMass=new TH1F("MissingMass_all","",missingmassparam);
	gHistoManager->Add(MissingMass,"Kinematics");
	{
		TH1F* hist=new TH1F("MissingMass_out_of_border","",missingmassparam);
		MissingMassDetailed.push_back(hist);
		gHistoManager->Add(hist,"Kinematics");
	}
	for(int i=1,N=P_Beam->GetNbinsX();i<=N;i++){
		int index=int(P_Beam->GetBinCenter(i)*1000);
		string name="MissingMass";
		name+=to_string(index);
		TH1F* hist=new TH1F(name.c_str(),"",missingmassparam);
		MissingMassDetailed.push_back(hist);
		gHistoManager->Add(hist,"Kinematics");
	}
#undef missingmassparam
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
	for(int i=0,n=ForwadrPlaneCount()-1;i<n;i++)
		EDepHist[i]->Fill(EDep(track,i+1),EDep(track,i));
	int stop_index=StoppingPlaneIndex(track);
	if(
		(
			(stop_index==ForwardPlaneIndex(kFTH1))||
			(stop_index==ForwardPlaneIndex(kFTH2))||
			(stop_index==ForwardPlaneIndex(kFTH3))||
			(stop_index==ForwardPlaneIndex(kFRH1))
		)
		&&UpperThresholdUpTo(track,stop_index)
		&&(EDep(track,stop_index)>0.01)
	){
		double Ek=0;
		if(ReconstructEkin(track,Ek)){
			for(int i=0,n=ForwadrPlaneCount()-1;i<n;i++)
				EDepFilteredHist[i]->Fill(EDep(track,i+1),EDep(track,i));
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
			TLorentzVector P_Total;{
				TLorentzVector P_Beam;
				TLorentzVector P_Target;
				P_Beam.SetVectM(p_beam,m_p);
				TVector3 ptarget;
				ptarget.SetMagThetaPhi(0,0,0);
				P_Target.SetVectM(ptarget,m_d);
				P_Total=P_Beam+P_Target;
			}
			TLorentzVector P_Missing=P_Total-P_He3;
			double missingmass=P_Missing.M();
			MissingMass->Fill(missingmass);
			DependenceOnPBeam->Fill(p_beam.Mag());
			{int index=0;
				for(int i=1,N=P_Beam->GetNbinsX();(i<=N)&&(index==0);i++){
					double min=P_Beam->GetBinLowEdge(i);
					if((p_beam.Mag()>=min)&&(p_beam.Mag()<(min+P_Beam->GetBinWidth(i))))
						index=i;
				}
				MissingMassDetailed[index]->Fill(missingmass);
			}
		}
	}
	return true;
}
bool He3eta::CentralTrackProcessing(WTrack* track,TVector3 &pbeam){
	return true;
}
void He3eta::EventPostProcessing(TVector3 &pbeam){}
