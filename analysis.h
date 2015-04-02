#ifndef QTBEVPHGEADUWGMY
#define QTBEVPHGEADUWGMY
#include <memory>
#include <functional>
#include <utility>
#include <vector>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <CAnalysisModule.hh>
#include <CDataManager.hh>
#include <FDEdep2Ekin.hh>
#include <CCardWDET.hh>
#include <Wasa.hh>
#include <CAnalysisModule.hh>
#include <REventWmcHeader.hh>
#include <REventHeader.hh>
#include <WTrackBank.hh>
#include <WVertexBank.hh>
#include <FDFTHTracks.hh>
#include <CDTracksSimple.hh>
enum TrackType{kFDN=1,kFDC=2,kCDN=11,kCDC=12};
enum ParticleType{
	kDummy=0,kGamma=1,kElectron=2,kPositron=3,kPi0=7,kPiPlus=8,kPiMinus=9,
	kNeutron=13,kProton=14,kEta=17,kDeuteron=45,kTriton=46,kHe3=49
};
class IAnalysis{
public:
	virtual ~IAnalysis();
	virtual void ProcessEvent()=0;
};
class Analysis:public virtual IAnalysis{
protected:
	Analysis();
public:
	virtual ~Analysis();
	virtual void ProcessEvent()override;
protected:
	virtual bool EventProcessingCondition()=0;
	virtual double PBeam()=0;
	virtual double EventWeight()=0;
	virtual void PrepareCheck();
	virtual void CheckParticleTrack(ParticleType type,double Ekin,double theta, double phi);

	virtual bool EventPreProcessing(TVector3 &pbeam)=0;
	virtual void EventPostProcessing(TVector3 &pbeam)=0;
	virtual bool TrackCountTrigger(int CinC,int NinC,int CinF)=0;
	virtual bool CentralFirst()=0;
	virtual bool ForwardTrackProcessing(WTrack* track,TVector3 &pbeam)=0;
	virtual bool CentralTrackProcessing(WTrack* track,TVector3 &pbeam)=0;

	WTrackBank *fTrackBankFD,*fTrackBankCD;
	CCardWDET *fDetectorTable;
	TH1F *P_Beam;
	vector<pair<ParticleType,double>> first_particles;
	vector<pair<ParticleType,double>> final_particles;
private:
	bool checkprepared;
};
template <class datatype,class reaction>
class CreateAnalysis:public virtual datatype,public virtual reaction{
public:
	CreateAnalysis():Analysis(),datatype(),reaction(){}
	virtual ~CreateAnalysis(){}
protected:
	virtual bool EventProcessingCondition()override{
		return datatype::EventProcessingCondition();
	}
	virtual double PBeam()override{
		return datatype::PBeam();
	}
	virtual double EventWeight()override{
		return datatype::EventWeight();
	}
	virtual void PrepareCheck()override{
		datatype::PrepareCheck();
	}
	virtual void CheckParticleTrack(ParticleType type,double Ekin,double theta, double phi)override{
		datatype::CheckParticleTrack(type,Ekin,theta,phi);
	}
	virtual bool EventPreProcessing(TVector3 &pbeam)override{
		return reaction::EventPreProcessing(pbeam);
	}
	virtual void EventPostProcessing(TVector3 &pbeam)override{
		reaction::EventPostProcessing(pbeam);
	}
	virtual bool TrackCountTrigger(int CinC,int NinC,int CinF)override{
		return reaction::TrackCountTrigger(CinC,NinC,CinF);
	}
	virtual bool CentralFirst()override{
		return reaction::CentralFirst();
	}
	virtual bool ForwardTrackProcessing(WTrack* track,TVector3 &pbeam)override{
		return reaction::ForwardTrackProcessing(track,pbeam);
	}
	virtual bool CentralTrackProcessing(WTrack* track,TVector3 &pbeam)override{
		return reaction::CentralTrackProcessing(track,pbeam);
	}
};
class MonteCarlo:public virtual Analysis{
public:
	MonteCarlo();
	virtual ~MonteCarlo();
protected:
	virtual bool EventProcessingCondition()override;
	virtual double PBeam()override;
	virtual double EventWeight()override;
	virtual void PrepareCheck()override;
	virtual void CheckParticleTrack(ParticleType type,double Ekin,double theta, double phi)override;
private:
	struct CheckHists{
	public:
		CheckHists(ParticleType t);
		ParticleType type;
		TH1F *Ekin, *Theta, *Phi;
	};
	vector<CheckHists> check;
	WTrackBank *fMCTrackBank;
	WVertexBank *fMCVertexBank;
	REventWmcHeader   *fEventHeader;
	
};
#endif
