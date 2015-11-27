// this file is distributed under 
// GPL v 3.0 license
#include <TCutG.h>
#include "../phys_constants.h"
#include "../kinematics.h"
#include "he3.h"
#include "detectors.h"
using namespace std;
He3_in_forward::He3_in_forward(double Q_lo,double Q_hi,unsigned int bins):Analysis(),ForwardDetectors(4),
	Reconstruction("Reconstruction",[this](){return 1000.0*Q_He3eta(PBeam());},bins,Q_lo,Q_hi),
	MissingMass("MissingMass",Reconstruction)
{
	AddLogSubprefix("He3");
	SubLog log=Log(LogDebug);
	AddParticleToFirstVertex(kHe3,m_3He);
	ForwardLayerCuts.push_back(
		TrackConditionSet("FRH1",Reconstruction).AddCondition("Stopping",[this](WTrack&track,vector<double>&){
			return (StopPlane(track)==kFRH1);
		}).AddCondition("He3Locus",[this](WTrack&track,vector<double>&){
			//Achtung - static
			static TCutG *cut=nullptr;
			if(cut==nullptr){
				cut=new TCutG("FRH1_cut",15);
				cut->SetVarX("FRH1");
				cut->SetVarY("FTH1");
				cut->SetPoint(15,0.019,0.025);
				cut->SetPoint(14,0.069,0.019);
				cut->SetPoint(13,0.121,0.016);
				cut->SetPoint(12,0.162,0.014);
				cut->SetPoint(11,0.206,0.013);
				cut->SetPoint(10,0.304,0.012);
				cut->SetPoint(9,0.4,0.012);
				cut->SetPoint(8,0.4,0.008);
				cut->SetPoint(7,0.298,0.008);
				cut->SetPoint(6,0.201,0.010);
				cut->SetPoint(5,0.141,0.012);
				cut->SetPoint(4,0.105,0.013);
				cut->SetPoint(3,0.061,0.017);
				cut->SetPoint(2,0.027,0.021);
				cut->SetPoint(1,0.019,0.025);
			}
			return cut->IsInside(EDep(track,kFRH1),EDep(track,kFTH1));
		}).AddParameter("E",[this](WTrack&track){
			static InterpolationBasedReconstruction energy("He3.E.FRH1"
				,[this](WTrack&track){return EDep(track,kFRH1);}
				,[this](WTrack&){return FromFirstVertex(kHe3).E;}
			);
			return energy.Reconstruct(track);
		})
	);
	ForwardLayerCuts.push_back(
		TrackConditionSet("FRH2",Reconstruction).AddCondition("Stopping",[this](WTrack&track,vector<double>&){
			return (StopPlane(track)==kFRH2);
		}).AddCondition("He3Locus",[this](WTrack&track,vector<double>&){
			vector<double> Ed={EDep(track,kFRH1),EDep(track,kFRH2)};
			return 
				(Ed[0]>(0.25-0.417*Ed[1]))&&(Ed[0]>(0.35-0.417*Ed[1]))
				&&(Ed[1]<0.22);
		}).AddParameter("E",[this](WTrack&track){
			//Achtung - static
			static InterpolationBasedReconstruction energy("He3.E.FRH2"
				,[this](WTrack&track){return EDep(track,kFRH1)+EDep(track,kFRH2);}
				,[this](WTrack&){return FromFirstVertex(kHe3).E;}
			);
			return energy.Reconstruct(track);
		})
	);
	Reconstruction.AddCondition("Theta_reconstruction_correct",[this](WTrack&track,vector<double>&P){
		//ToDo: replace by more reasonable condition
		return 0.125!=track.Theta();//Magic number taken from framework
	}).AddCondition("Edep_cuts",[this](WTrack&track,vector<double>&P){
		ForwardDetectorTrackMarker(0,track);
		bool res=false;
		for(auto&layer:ForwardLayerCuts)res|=layer.Check(track,P);
		return res;
	}).AddParameter("Theta",[this](WTrack&track){
		ForwardDetectorTrackMarker(1,track);
		auto Th_m=[this](WTrack&track){return track.Theta();};
		auto Th_t=[this](WTrack&){return FromFirstVertex(kHe3).Th;};
		//Achtung - static
		static InterpolationBasedReconstruction He3_theta("He3.th",Th_m,Th_t);
		return He3_theta.Reconstruct(track);
	}).AddParameter("Phi",[this](WTrack&track){
		auto Ph_m=[this](WTrack&track){return NormPhi(track.Phi());};
		auto Ph_t=[this](WTrack&){return FromFirstVertex(kHe3).Phi;};
		//Achtung - static
		static InterpolationBasedReconstruction He3_phi("He3.phi",Ph_m,Ph_t);
		return He3_phi.Reconstruct(track);
	}).AddCondition("Reconstructed",[this](WTrack&track,vector<double>&P){
		return isfinite(P[0])&&isfinite(P[1])&&isfinite(P[2]);
	}).AddCondition("Additional",[this](WTrack&track,vector<double>&P){
		ForwardDetectorTrackMarker(2,track);
		for(auto condition:this->AdditionalConditions)
			if(condition(track)==false)
				return false;
		ForwardDetectorTrackMarker(3,track);
		return true;
	});
	MissingMass.Setup([this](vector<double>&He3){
		double p=sqrt(He3[0]*(He3[0]+2*m_3He));
		TVector3 p_He3;
		p_He3.SetMagThetaPhi(p,He3[1],He3[2]);
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
		return P_Missing.M();
	},300,m_pi0-0.05,m_eta+0.05);
	AddTrackProcessing(kFDC,[this](WTrack&track){
		vector<double> He3params;
		if(Reconstruction.Check(track,He3params))
			MissingMass.AcceptEvent(He3params);
	});
}
He3_in_forward::~He3_in_forward(){}
void He3_in_forward::AddCondition(ConditionTrackDependent condition){
	AdditionalConditions.push_back(condition);
}
bool He3_in_forward::EventPreProcessing(){
	Reconstruction.ReferenceEvent();
	return true;
}
bool He3_in_forward::TrackCountTrigger(int CinC,int NinC,int CinF,int NinF){
	return CinF>0;
}
void He3_in_forward::EventPostProcessing(){}

He3_Modification_for_eta::He3_Modification_for_eta():He3_in_forward(0.0,30.0,12){
	AddCondition([this](WTrack&track){
		if(StopPlane(track)!=kFRH1)
			return false;
		double E=EDep(track,kFRH1);
		return (E>0.08)&&(E<0.22);
	});
}
MC_He3eta::MC_He3eta(){
	AddParticleToFirstVertex(kEta,m_eta);
}
MC_He3pi0::MC_He3pi0(){
	AddParticleToFirstVertex(kPi0,m_pi0);
}
