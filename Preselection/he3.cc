// this file is distributed under 
// GPL v 3.0 license
#include <TCutG.h>
#include "he3.h"
#include "detectors.h"
#include "../General/phys_constants.h"
#include "geometry.h"
using namespace std;
He3_in_forward::He3_in_forward():Analysis(),ForwardDetectors(2),
	Cut("Cuts",[this](){return PBeam();},20,p_he3_eta_threshold,p_beam_hi),
	He3_Ekin("He3.E.FTH1nFRH1",
		 [this](WTrack&track){return EDep(track,kFRH1)+EDep(track,kFTH1);},
		 [this](WTrack&){return FromFirstVertex(kHe3).E;}
	),
	He3_theta("He3.th",
		  [this](WTrack&track){return track.Theta();},
		  [this](WTrack&){return FromFirstVertex(kHe3).Th;}
	),
	He3_phi("He3.phi",
		[this](WTrack&track){return NormPhi(track.Phi());},
		[this](WTrack&){return FromFirstVertex(kHe3).Phi;}
	),
	MissingMass("MissingMass",Cut){
	AddLogSubprefix("He3");
	SubLog log=Log(LogDebug);
	AddParticleToFirstVertex(kHe3,m_3He);
	AddParticleToFirstVertex(kEta,m_eta);
	Cut.AddCondition("StoppingPlane",[this](WTrack&track,vector<double>&){
		return (StopPlane(track)==kFRH1);
	}).AddCondition("SP2cut",[this](WTrack&track,vector<double>&){
		using namespace PlaneGeometry;
		static Polygon<double> mycut;
		if(mycut.size()==0){
			mycut
				<<make_pair(0.000,0.033)
				<<make_pair(0.068,0.023)
				<<make_pair(0.146,0.017)
				<<make_pair(0.223,0.016)
				<<make_pair(0.284,0.013)
				<<make_pair(0.214,0.009)
				<<make_pair(0.139,0.009)
				<<make_pair(0.093,0.013)
				<<make_pair(0.047,0.016)
				<<make_pair(0.000,0.024)
				<<make_pair(0.000,0.033);
		}
		return mycut(make_pair(EDep(track,kFRH1),EDep(track,kFTH1)));
	}).AddCondition("ThetaCut",[this](WTrack&track,vector<double>&){
		return ((track.Theta()<0.1245)||(track.Theta()>0.1255));
	}).AddParameter("E",[this](WTrack&track){
		return He3_Ekin.Reconstruct(track);
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
}
He3_in_forward::~He3_in_forward(){}
bool He3_in_forward::EventPreProcessing(){
	Cut.ReferenceEvent();
	return true;
}
bool He3_in_forward::TrackCountTrigger(int CinC,int NinC,int CinF){return CinF>0;}
bool He3_in_forward::CentralFirst(){return false;}
bool He3_in_forward::ForwardTrackProcessing(WTrack&track){
	SubLog log=Log(LogDebug);
	ForwardDetectorTrackMarker(0,track);
	vector<double> He3;
	if(Cut.Check(track,He3)){
		ForwardDetectorTrackMarker(1,track);
		MissingMass.AcceptEvent(He3);
	}
	return true;
}
bool He3_in_forward::CentralTrackProcessing(WTrack&){return true;}
void He3_in_forward::EventPostProcessing(){}
