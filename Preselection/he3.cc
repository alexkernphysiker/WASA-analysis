// this file is distributed under 
// GPL v 3.0 license
#include <TCutG.h>
#include "he3.h"
#include "detectors.h"
using namespace std;
He3_in_forward::He3_in_forward():Analysis(),ForwardDetectors(2),
	Cut("Cuts",[this](){return PBeam();},20,p_he3_eta_threshold,p_beam_hi),
	CutFTH1("FTH1",Cut),CutFRH1("FRH1",Cut),
	MissingMass("MissingMass",Cut)
{
	AddLogSubprefix("He3");
	SubLog log=Log(LogDebug);
	AddParticleToFirstVertex(kHe3,m_3He);

	chargehist=new TH1F("Charge","",16,0,15);
	gHistoManager->Add(chargehist,"Especial");
	
	Th_m=[this](WTrack&track){return track.Theta();};
	Th_t=[this](WTrack&){return FromFirstVertex(kHe3).Th;};
	Ph_m=[this](WTrack&track){return NormPhi(track.Phi());};
	Ph_t=[this](WTrack&){return FromFirstVertex(kHe3).Phi;};
	Elow_m=[this](WTrack&track){return EDep(track,kFTH1);};
	Ehi_m=[this](WTrack&track){return EDep(track,kFRH1)+EDep(track,kFTH1);};
	E_t=[this](WTrack&){return FromFirstVertex(kHe3).E;};

	He3_theta=InterpolationBasedReconstruction("He3.th",Th_m,Th_t);
	He3_phi=InterpolationBasedReconstruction("He3.phi",Ph_m,Ph_t);
	He3_Ekin.push_back(InterpolationBasedReconstruction("He3.E.FTH1",Elow_m,E_t));
	He3_Ekin.push_back(InterpolationBasedReconstruction("He3.E.FRH1",Ehi_m,E_t));
	CutFRH1.AddCondition("StoppingFRH1",[this](WTrack&track,vector<double>&){
		return (StopPlane(track)==kFRH1);
	}).AddCondition("SP2cut",[this](WTrack&track,vector<double>&){
		static TCutG *cut=nullptr;
		if(cut==nullptr){
			cut=new TCutG("FRH1_cut",12);
			cut->SetVarX("FRH1");
			cut->SetVarY("FTH1");
			cut->SetPoint(11,0.000,0.033);
			cut->SetPoint(10,0.068,0.023);
			cut->SetPoint(9,0.146,0.017);
			cut->SetPoint(8,0.223,0.016);
			cut->SetPoint(7,0.284,0.013);
			cut->SetPoint(6,0.214,0.009);
			cut->SetPoint(5,0.139,0.009);
			cut->SetPoint(4,0.093,0.013);
			cut->SetPoint(3,0.047,0.016);
			cut->SetPoint(2,0.000,0.024);
			cut->SetPoint(1,0.000,0.033);
		}
		return cut->IsInside(EDep(track,kFRH1),EDep(track,kFTH1));
	}).AddParameter("E",[this](WTrack&track){
		return He3_Ekin[1].Reconstruct(track);
	});
	CutFTH1.AddCondition("StoppingFTH1",[this](WTrack&track,vector<double>&){
		return (StopPlane(track)==kFTH1);
	}).AddCondition("FWC2cut",[this](WTrack&track,vector<double>&){
		return EDep(track,kFWC2)>=0.011;
	}).AddParameter("E",[this](WTrack&track){
		return He3_Ekin[0].Reconstruct(track);
	});
	Cut.AddCondition("ReconstructionCondition",[this](WTrack&track,vector<double>&){
		auto type=track.ParticleType();
		chargehist->Fill(type);
		bool passes=(type==kFDC);
		if(passes)debug_yes(track);
		else debug_no(track);
		return passes;
	}).AddCondition("Edep_cuts",[this](WTrack&track,vector<double>&P){
		auto one=CutFRH1.Check(track,P),two=CutFTH1.Check(track,P);
		return one||two;//for we could see correct numbers of events on histograms
	}).AddParameter("Theta",[this](WTrack&track){
		return He3_theta.Reconstruct(track);
	}).AddParameter("Phi",[this](WTrack&track){
		return He3_phi.Reconstruct(track);
	}).AddCondition("Reconstructed",[this](WTrack&,vector<double>&P){
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
	},400,0.4,0.6);
	AddTrackProcessing(make_pair(kFDC,[this](WTrack&track){
		SubLog log=Log(LogDebug);
		ForwardDetectorTrackMarker(0,track);
		vector<double> He3;
		if(Cut.Check(track,He3)){
			ForwardDetectorTrackMarker(1,track);
			MissingMass.AcceptEvent(He3);
		}
	}));
}
He3_in_forward::~He3_in_forward(){}
bool He3_in_forward::EventPreProcessing(){
	Cut.ReferenceEvent();
	return true;
}
bool He3_in_forward::TrackCountTrigger(int CinC,int NinC,int CinF,int NinF){
	return CinF>0;
}
void He3_in_forward::EventPostProcessing(){}
void He3_in_forward::debug_yes(WTrack&){}
void He3_in_forward::debug_no(WTrack&){}

He3_mc_debug::He3_mc_debug():He3_in_forward(),Yes("ReconstructionConditionYes"),No("ReconstructionConditionNo"){
	auto prepare=[this](Debug2DSpectraSet&D){
		Axis Ekin={.value=E_t,.from=0,.to=0.5,.bins=500};
		Axis Ehi={.value=Ehi_m,.from=0,.to=0.3,.bins=300};
		Axis Tht={.value=Th_t,.from=0,.to=0.3,.bins=300};
		Axis Thm={.value=Th_m,.from=0.09,.to=0.16,.bins=70};
		Axis Pht={.value=Ph_t,.from=0,.to=7,.bins=700};
		Axis Phm={.value=Ph_m,.from=0,.to=7,.bins=700};
		D.Add("Vertex_E_vs_Theta",Ekin,Tht);
		D.Add("Vertex_E_vs_Phi",Ekin,Pht);
		D.Add("Vertex_Theta_vs_Phi",Tht,Pht);
		D.Add("Detector_E_vs_Th",Ehi,Thm);
		D.Add("Detector_E_vs_Ph",Ehi,Phm);
		D.Add("Detector_Th_vs_Ph",Thm,Phm);
	};
	prepare(Yes);
	prepare(No);
}
He3_mc_debug::~He3_mc_debug(){}
void He3_mc_debug::debug_yes(WTrack&track){
	Yes.CatchState(track);
}
void He3_mc_debug::debug_no(WTrack&track){
	No.CatchState(track);
}

He3eta::He3eta():He3_mc_debug(){AddParticleToFirstVertex(kEta,m_eta);}
He3pi0::He3pi0():He3_mc_debug(){AddParticleToFirstVertex(kPi0,m_pi0);}
