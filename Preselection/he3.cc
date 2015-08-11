// this file is distributed under 
// GPL v 3.0 license
#include "he3.h"
#include "detectors.h"
#include "../General/phys_constants.h"
using namespace std;
He3_at_FRH1::He3_at_FRH1():Analysis(),
He3_Ekin("He3.E",
	[this](WTrack&&track){
		double res=0;
		for(int i=0,n=ForwadrPlaneCount();i<(n-1);i++)
			res+=EDep(static_right(track),i);
		return res;
	},
	[this](WTrack&&track){return FromFirstVertex(kHe3).E;}
),
He3_theta("He3.th",
	[this](WTrack&&track){return track.Theta();},
	[this](WTrack&&track){return FromFirstVertex(kHe3).Th;}
),
He3_phi("He3.phi",
	[this](WTrack&&track){return NormPhi(track.Phi());},
	[this](WTrack&&track){return FromFirstVertex(kHe3).Phi;}
){
	AddLogSubprefix("He3");
	SubLog log=Log(LogDebug);
	AddParticleToFirstVertex(kHe3,m_3He);
	AddParticleToFirstVertex(kEta,m_eta);
	for(int i=0,n=ForwadrPlaneCount()-1;i<n;i++){
		string histname=ForwardPlaneName(i)+"_vs_"+ForwardPlaneName(i+1);
		TH2F* hist;
		auto makehist=[this,i,&hist,&histname,&log](vector<TH2F*> &vect){
			hist=new TH2F(histname.c_str(),"",128,0,UpperByIndex(i+1),128,0,UpperByIndex(i));
			vect.push_back(hist);
			log<<"Adding Edep 2D-histogram";
			log<<histname;
			gHistoManager->Add(hist,"E_dep");
		};
		makehist(EDepHist);
		histname=histname+"_filtered";
		makehist(EDepFilteredHist);
	}
	P_Beam=new TH1F("AllEventsOnPBeam","",20,p_he3_eta_threshold,p_beam_hi);
	gHistoManager->Add(P_Beam,"EventsCount");
	DependenceOnPBeam=new TH1F("FilteredEventsOnPBeam","",20,p_he3_eta_threshold,p_beam_hi);
	gHistoManager->Add(DependenceOnPBeam,"EventsCount");
	#define missingmassparam 400,0.4,0.6
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
		log<<"Add missing mass histogram";
		gHistoManager->Add(hist,"Kinematics");
	}
	#undef missingmassparam
}
He3_at_FRH1::~He3_at_FRH1(){}
bool He3_at_FRH1::EventPreProcessing(){
	P_Beam->Fill(PBeam());
	return true;
}
bool He3_at_FRH1::TrackCountTrigger(int CinC,int NinC,int CinF){return CinF>0;}
bool He3_at_FRH1::CentralFirst(){return false;}
bool He3_at_FRH1::ForwardTrackProcessing(WTrack&& track){
	SubLog log=Log(LogDebug);
	for(int i=0,n=ForwadrPlaneCount()-1;i<n;i++)
		EDepHist[i]->Fill(EDep(static_right(track),i+1),EDep(static_right(track),i));
	if(
		(StopPlane(static_right(track))==kFRH1)
		&&(EDep(static_right(track),kFWC1)>0.0055)
		&&(EDep(static_right(track),kFWC1)<0.015)
		&&(EDep(static_right(track),kFWC2)>0.0055)
		&&(EDep(static_right(track),kFWC2)<0.015)
		&&(EDep(static_right(track),kFTH1)>0.009)
		&&(EDep(static_right(track),kFTH1)<0.030)
		&&(EDep(static_right(track),kFRH1)>0.05)
		&&(EDep(static_right(track),kFRH1)<0.30)
		&&((track.Theta()<0.1245)||(track.Theta()>0.1255))
	){
		log<<"stopping plane index condition passed";
		for(int i=0,n=ForwadrPlaneCount()-1;i<n;i++)
			EDepFilteredHist[i]->Fill(EDep(static_right(track),i+1),EDep(static_right(track),i));
		double Ek,th,phi;
		bool c1=He3_Ekin.Reconstruct(Ek,static_right(track)),
			c2=He3_theta.Reconstruct(th,static_right(track)),
			c3=He3_phi.Reconstruct(phi,static_right(track));
		if(c1&&c2&&c3){
			log<<"reconstruction successful";
			double p=sqrt(Ek*(Ek+2*m_3He));
			TVector3 p_He3;
			p_He3.SetMagThetaPhi(p,th,phi);
			TLorentzVector P_He3;
			P_He3.SetVectM(p_He3,m_3He);
			TLorentzVector P_Total;{
				TVector3 p_beam;
				p_beam.SetMagThetaPhi(PBeam(),0,0);
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
			DependenceOnPBeam->Fill(PBeam());
			{int index=0;
				for(int i=1,N=P_Beam->GetNbinsX();(i<=N)&&(index==0);i++){
					double min=P_Beam->GetBinLowEdge(i);
					if((PBeam()>=min)&&(PBeam()<(min+P_Beam->GetBinWidth(i))))
						index=i;
				}
				MissingMassDetailed[index]->Fill(missingmass);
			}
		}
	}
	return true;
}
bool He3_at_FRH1::CentralTrackProcessing(WTrack&& track){return true;}
void He3_at_FRH1::EventPostProcessing(){}
