#ifndef QTBEVPHGEADUWGMY
#define QTBEVPHGEADUWGMY
#include <memory>
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
	virtual void ProcessEvent()final;
protected:
	virtual bool EventProcessingCondition()=0;
	virtual double PBeam()=0;
	virtual double EventWeight()=0;
	virtual void PrepareCheck()=0;
	virtual void CheckParticleTrack(ParticleType type,double Ekin,double theta, double phi)=0;

	virtual bool EventPreProcessing(TVector3 &pbeam)=0;
	virtual void EventPostProcessing(TVector3 &pbeam)=0;
	virtual bool TrackCountTrigger(int CinC,int NinC,int CinF)=0;
	virtual bool CentralFirst()=0;
	virtual bool ForwardTrackProcessing(WTrack* track,TVector3 &pbeam)=0;
	virtual bool CentralTrackProcessing(WTrack* track,TVector3 &pbeam)=0;

	FDFTHTracks* TrackFinderFD;
	CDTracksSimple* CDTrackFinder;
	WTrackBank *fTrackBankFD,*fTrackBankCD;
	CCardWDET *fDetectorTable;
	vector<pair<ParticleType,double>> first_particles;
	vector<pair<ParticleType,double>> final_particles;
private:
	bool checkprepared;
};
template <class datatype,class reaction>
class CustomAnalysis:public virtual datatype,public virtual reaction{
public:
	CustomAnalysis():Analysis(),datatype(),reaction(){}
	virtual ~CustomAnalysis(){}
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
#endif
