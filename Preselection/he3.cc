// this file is distributed under 
// GPL v 3.0 license
#include <TCutG.h>
#include "he3.h"
#include "detectors.h"
#include "../General/phys_constants.h"
using namespace std;
He3_in_forward::He3_in_forward():Analysis(),ForwardDetectors(2),
	Cut("Cuts",[this](){return PBeam();},20,p_he3_eta_threshold,p_beam_hi),
	He3_Ekin("He3.E.FTH1nFRH1",
		 [this](const WTrack&track){return EDep(track,kFRH1)+EDep(track,kFTH1);},
		 [this](const WTrack&){return FromFirstVertex(kHe3).E;}
	),
	He3_theta("He3.th",
		  [this](const WTrack&track){return track.Theta();},
		  [this](const WTrack&){return FromFirstVertex(kHe3).Th;}
	),
	He3_phi("He3.phi",
		[this](const WTrack&track){return NormPhi(track.Phi());},
		[this](const WTrack&){return FromFirstVertex(kHe3).Phi;}
	),
	MissingMass("MissingMass",Cut){
	AddLogSubprefix("He3");
	SubLog log=Log(LogDebug);
	AddParticleToFirstVertex(kHe3,m_3He);
	AddParticleToFirstVertex(kEta,m_eta);
	Cut.AddCondition("StoppingPlane",[this](const WTrack&track,vector<double>&){
		return (StopPlane(track)==kFRH1);
	}).AddCondition("SP2cut",[this](const WTrack&track,vector<double>&){
		static TCutG *cut=nullptr;
		if(cut==nullptr){
			cut=new TCutG("FRH1_cut",11);
			cut->SetVarX("FRH1");
			cut->SetVarY("FTH1");
			cut->SetPoint(0,0.000,0.033);
			cut->SetPoint(2,0.068,0.023);
			cut->SetPoint(3,0.146,0.017);
			cut->SetPoint(4,0.223,0.016);
			cut->SetPoint(5,0.284,0.013);
			cut->SetPoint(6,0.214,0.009);
			cut->SetPoint(7,0.139,0.009);
			cut->SetPoint(8,0.093,0.013);
			cut->SetPoint(9,0.047,0.016);
			cut->SetPoint(10,0.000,0.024);
		}
		return cut->IsInside(EDep(track,kFRH1),EDep(track,kFTH1));
	}).AddCondition("ThetaCut",[this](const WTrack&track,vector<double>&){
		return ((track.Theta()<0.1245)||(track.Theta()>0.1255));
	}).AddParameter("E",[this](const WTrack&track){
		return He3_Ekin.Reconstruct(track);
	}).AddParameter("Theta",[this](const WTrack&track){
		return He3_theta.Reconstruct(track);
	}).AddParameter("Phi",[this](const WTrack&track){
		return He3_phi.Reconstruct(track);
	}).AddCondition("Reconstructed",[this](const WTrack&,vector<double>&P){
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
}
He3_in_forward::~He3_in_forward(){}
bool He3_in_forward::EventPreProcessing(){
	Cut.ReferenceEvent();
	return true;
}
bool He3_in_forward::TrackCountTrigger(int CinC,int NinC,int CinF){return CinF>0;}
bool He3_in_forward::CentralFirst(){return false;}
bool He3_in_forward::ForwardTrackProcessing(const WTrack&track){
	SubLog log=Log(LogDebug);
	ForwardDetectorTrackMarker(0,track);
	vector<double> He3;
	if(Cut.Check(track,He3)){
		ForwardDetectorTrackMarker(1,track);
		MissingMass.AcceptEvent(He3);
	}
	return true;
}
bool He3_in_forward::CentralTrackProcessing(const WTrack&){return true;}
void He3_in_forward::EventPostProcessing(){}
