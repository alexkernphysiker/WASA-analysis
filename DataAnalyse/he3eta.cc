#include "he3eta.h"
#include "detectors.h"
#include "../General/phys_constants.h"
using namespace std;
He3eta::He3eta():Analysis(),ForwardDetectorRoutines("3He"){
	AddSubprefix("He3eta");
	SubLog log=getSubLog("Constructor");
	first_particles.push_back(make_pair(kHe3,m_3He));
	first_particles.push_back(make_pair(kEta,m_eta));
	final_particles.push_back(make_pair(kHe3,m_3He));
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
	SetGettableFunction([](ForwardDetectorPlane StopPlane,int& Edep2Ekin_table){
		if(StopPlane==kFVH || StopPlane==kFTH3 || StopPlane==kFRH3)
			return false;
		Edep2Ekin_table=int(StopPlane);
		return true;
	});
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
He3eta::~He3eta(){
	SubLog log=getSubLog("Destructor");
}
bool He3eta::EventPreProcessing(TVector3 &pbeam){
	P_Beam->Fill(pbeam.Mag());
	return true;
}
bool He3eta::TrackCountTrigger(int CinC,int NinC,int CinF){
	return true;
}
bool He3eta::CentralFirst(){
	return false;
}
bool He3eta::ForwardTrackProcessing(WTrack* track,TVector3 &p_beam){
	SubLog log=getSubLog("ForwardTrackProcessing");
	for(int i=0,n=ForwadrPlaneCount()-1;i<n;i++)
		EDepHist[i]->Fill(EDep(track,i+1),EDep(track,i));
	if(
		(StopPlane(track)==kFRH1)
		&&(EDep(track,kFTH1)>0.01)
	){
		double Ek=0;
		log<<"stopping plane index condition passed";
		if(ReconstructEkin(track,Ek)){
			log<<"reconstruction successful";
			for(int i=0,n=ForwadrPlaneCount()-1;i<n;i++)
				EDepFilteredHist[i]->Fill(EDep(track,i+1),EDep(track,i));
			double p=sqrt(Ek*(Ek+2*m_3He));
			double theta=track->Theta();
			double phi=track->Phi();
			double magnetic_field=10;//kG
			//ToDo: obtain magnetic field
			phi+= TrackFinderFD->GetPhiCorrection(p,theta,magnetic_field,c_He);
			CheckParticleTrack(kHe3,Ek,theta,phi);
			TVector3 p_He3;
			p_He3.SetMagThetaPhi(p,theta,phi);
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
		}log.Message(LogWarning,"Could not reconstruct energy");
	}
	return true;
}
bool He3eta::CentralTrackProcessing(WTrack* track,TVector3 &pbeam){
	return true;
}
void He3eta::EventPostProcessing(TVector3 &pbeam){}
