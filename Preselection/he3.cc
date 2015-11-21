// this file is distributed under 
// GPL v 3.0 license
#include <TCutG.h>
#include "../phys_constants.h"
#include "he3.h"
#include "detectors.h"
using namespace std;
He3_in_forward::He3_in_forward():Analysis(),ForwardDetectors(2),
	Reconstruction("Reconstruction",[this](){return PBeam();},beam_momenta_bins,p_beam_low,p_beam_hi),
	MissingMass("MissingMass",Reconstruction)
{
	AddLogSubprefix("He3");
	SubLog log=Log(LogDebug);
	AddParticleToFirstVertex(kHe3,m_3He);

	ForwardLayerCuts.push_back(
		TrackConditionSet("FTH1",Reconstruction).AddCondition("Stopping",[this](WTrack&track,vector<double>&){
			return (StopPlane(track)==kFTH1);
		}).AddCondition("FWC2cut",[this](WTrack&track,vector<double>&){
			return EDep(track,kFWC2)>=0.011;
		}).AddParameter("E",[this](WTrack&track){
			static InterpolationBasedReconstruction energy("He3.E.FTH1"
				,[this](WTrack&track){return EDep(track,kFTH1);}
				,[this](WTrack&){return FromFirstVertex(kHe3).E;}
			);
			return energy.Reconstruct(track);
		})
	);
	ForwardLayerCuts.push_back(
		TrackConditionSet("FRH1",Reconstruction).AddCondition("Stopping",[this](WTrack&track,vector<double>&){
			return (StopPlane(track)==kFRH1);
		}).AddCondition("SP2cut",[this](WTrack&track,vector<double>&){
			//Achtung - static
			static TCutG *cut=nullptr;
			if(cut==nullptr){
				cut=new TCutG("FRH1_cut",13);
				cut->SetVarX("FRH1");
				cut->SetVarY("FTH1");
				cut->SetPoint(12,0.000,0.033);
				cut->SetPoint(11,0.068,0.023);
				cut->SetPoint(10,0.146,0.017);
				cut->SetPoint(9,0.223,0.016);
				cut->SetPoint(8,0.4,0.012);
				cut->SetPoint(7,0.4,0.008);
				cut->SetPoint(6,0.233,0.008);
				cut->SetPoint(5,0.139,0.009);
				cut->SetPoint(4,0.093,0.013);
				cut->SetPoint(3,0.047,0.016);
				cut->SetPoint(2,0.000,0.024);
				cut->SetPoint(1,0.000,0.033);
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
		}).AddCondition("Linear_cuts",[this](WTrack&track,vector<double>&){
			// two linear cuts
			return
				(EDep(track,kFRH1)>(0.25-0.417*EDep(track,kFRH2)))
				&&(EDep(track,kFRH1)<(0.35-0.417*EDep(track,kFRH2)))
				&&(EDep(track,kFRH2)<0.22);
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
		return 0.125!=track.Theta();//Magic number taken from framework
	}).AddCondition("Edep_cuts",[this](WTrack&track,vector<double>&P){
		ForwardDetectorTrackMarker(0,track);
		bool res=false;
		for(auto&layer:ForwardLayerCuts)res|=layer.Check(track,P);
		return res;
	}).AddParameter("Theta",[this](WTrack&track){
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
		ForwardDetectorTrackMarker(1,track);
		return isfinite(P[0])&&isfinite(P[1])&&isfinite(P[2]);
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
	},600,m_pi0-0.05,m_eta+0.05);
	AddTrackProcessing(kFDC,[this](WTrack&track){
		vector<double> He3params;
		if(Reconstruction.Check(track,He3params))
			MissingMass.AcceptEvent(He3params);
	});
}
He3_in_forward::~He3_in_forward(){}
bool He3_in_forward::EventPreProcessing(){
	Reconstruction.ReferenceEvent();
	return true;
}
bool He3_in_forward::TrackCountTrigger(int CinC,int NinC,int CinF,int NinF){
	return CinF>0;
}
void He3_in_forward::EventPostProcessing(){}

He3eta::He3eta():He3_in_forward(){AddParticleToFirstVertex(kEta,m_eta);}
He3pi0::He3pi0():He3_in_forward(){AddParticleToFirstVertex(kPi0,m_pi0);}
