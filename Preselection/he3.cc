// this file is distributed under 
// GPL v 3.0 license
#include "he3.h"
#include "detectors.h"
#include "../General/phys_constants.h"
using namespace std;
He3_at_FRH1::He3_at_FRH1():Analysis(),ForwardDetectors(2),
He3_Ekin("He3.E.FTH1nFRH1",
	[this](WTrack&&track){
		return EDep(static_right(track),kFRH1)+EDep(static_right(track),kFTH1);
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
),
Cuts("Events",[this](){return PBeam();},20,p_he3_eta_threshold,p_beam_hi),
MissingMass("MissingMass",static_right(Cuts),
	[this](vector<double>&He3){
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
	},400,0.4,0.6)
{
	AddLogSubprefix("He3");
	SubLog log=Log(LogDebug);
	AddParticleToFirstVertex(kHe3,m_3He);
	AddParticleToFirstVertex(kEta,m_eta);
	Cuts.AddCondition("StoppingPlane",[this](WTrack&&track,vector<double>&){
		return (StopPlane(static_right(track))==kFRH1);
	}).AddCondition("SP2cut",[this](WTrack&&track,vector<double>&){
		return (EDep(static_right(track),kFTH1)>0.009)&&(EDep(static_right(track),kFTH1)<0.030)
			&&(EDep(static_right(track),kFRH1)>0.05)&&(EDep(static_right(track),kFRH1)<0.30);
	}).AddCondition("ThetaCut",[this](WTrack&&track,vector<double>&){
		return ((track.Theta()<0.1245)||(track.Theta()>0.1255));
	}).AddParameter("E",[this](WTrack&&track){
		return He3_Ekin.Reconstruct(static_right(track));
	}).AddParameter("Theta",[this](WTrack&&track){
		return He3_theta.Reconstruct(static_right(track));
	}).AddParameter("Phi",[this](WTrack&&track){
		return He3_phi.Reconstruct(static_right(track));
	}).AddCondition("Reconstructed",[this](WTrack,vector<double>&P){
		return isfinite(P[0])&&isfinite(P[1])&&isfinite(P[2]);
	});
}
He3_at_FRH1::~He3_at_FRH1(){}
bool He3_at_FRH1::EventPreProcessing(){
	Cuts.ReferenceEvent();
	return true;
}
bool He3_at_FRH1::TrackCountTrigger(int CinC,int NinC,int CinF){return CinF>0;}
bool He3_at_FRH1::CentralFirst(){return false;}
bool He3_at_FRH1::ForwardTrackProcessing(WTrack&& track){
	SubLog log=Log(LogDebug);
	ForwardDetectorTrackMarker(0,static_right(track));
	vector<double> He3;
	if(Cuts.Check(static_right(track),He3)){
		ForwardDetectorTrackMarker(1,static_right(track));
		MissingMass.AcceptEvent(He3);
	}
	return true;
}
bool He3_at_FRH1::CentralTrackProcessing(WTrack&& track){return true;}
void He3_at_FRH1::EventPostProcessing(){}
