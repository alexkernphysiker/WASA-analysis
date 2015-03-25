#include "reactions.h"
He3eta_gg::He3eta_gg():Analysis(){
	MC_beam_momenta.push_back(make_pair(kHe3,m_3He));
	MC_beam_momenta.push_back(make_pair(kEta,m_eta));
	FRH1vsFRH2=new TH2F("FRH1vsFRH2","",256,0,0.5,256,0,0.5);
	gHistoManager->Add(FRH1vsFRH2,"FRH1vsFRH2");
}
He3eta_gg::~He3eta_gg(){}
bool He3eta_gg::EventPreProcessing(TVector3 &pbeam){
	return true;
}
bool He3eta_gg::TrackCountTrigger(int CinC,int NinC,int CinF){
	return true;
}
bool He3eta_gg::CentralFirst(){return false;}
bool He3eta_gg::ForwardTrackProcessing(WTrack* track,TVector3 &pbeam){
	const int fwd_count=5; 
	double fwd_thresholds[]={0.004,0.0025,0.0025,0.0035,0.004};
	ForwardDetectorPlane fwd_planes[]={kFRH1,kFRH2,kFRH3,kFRH4,kFRH5};
	auto thresholds=[track](){
		bool res=true;
		double thresholds[]={0.00018,0.00018,0.0015,0.00032,0.0003};
		ForwardDetectorPlane thr_planes[]={kFWC1,kFWC2,kFTH1,kFTH2,kFTH3};
		for(int i=0;i<5;i++)
			track->Edep(thr_planes[i])>thresholds[i];
		return res;
	};
	auto forward_stop_plane=[track,fwd_planes,fwd_thresholds](unsigned int stopindex){
		bool res=true;
		for(int i=0;i<=stopindex;i++)
			res&=track->Edep(fwd_planes[i])>fwd_thresholds[i];
		for(int i=stopindex+1;i<fwd_count;i++)
			res&=track->Edep(fwd_planes[i])<=fwd_thresholds[i];
		return res;
	};
	bool stop_condition=false;
	for(int i=0;i<fwd_count;i++)
		stop_condition|=forward_stop_plane(i);
	if(thresholds()&&stop_condition){
		double EdepFRHtot=0;
		for(int i=0;i<fwd_count;i++)
			EdepFRHtot+=track->Edep(fwd_planes[i]);
		//ToDo: preselection condition for Edep
		TLorentzVector P_b;
		P_b.SetVectM(pbeam,m_p);
		TVector3 vec_target;
		vec_target.SetMagThetaPhi(0,0,0);
		TLorentzVector P_t;
		P_t.SetVectM(vec_target,m_d);
		TLorentzVector P_tot=P_b+P_t;
	}
	return true;
}
bool He3eta_gg::CentralTrackProcessing(WTrack* track,TVector3 &pbeam){return true;}
void He3eta_gg::EventPostProcessing(TVector3 &pbeam){}