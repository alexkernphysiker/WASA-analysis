// this file is distributed under 
// GPL v 3.0 license
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
#include "log.h"
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
class Analysis:public virtual IAnalysis,public virtual Logger{
protected:
	Analysis();
public:
	virtual ~Analysis();
	virtual void ProcessEvent()final;
protected:
	virtual void PrepairForEventAnalysis()=0;
	virtual bool EventProcessingCondition()=0;

	virtual bool EventPreProcessing()=0;
	virtual void EventPostProcessing()=0;
	virtual bool TrackCountTrigger(int CinC,int NinC,int CinF,int NinF)=0;

	FDFTHTracks* TrackFinderFD;
	CDTracksSimple* CDTrackFinder;
	WTrackBank *fTrackBankFD,*fTrackBankCD;
	CCardWDET *fDetectorTable;
	double PBeam();
	void CachePBeam(double value);
	void AddParticleToFirstVertex(ParticleType type,double mass);
	struct Kinematic{
		Kinematic();
		double E,Th,Phi;
	};
	typedef pair<TrackType,function<void(WTrack&)>> TrackProcessing;
	void AddTrackProcessing(TrackProcessing&&proc);
	Kinematic&FromFirstVertex(ParticleType type);
	void ForFirstVertex(std::function<void(ParticleType,double,Kinematic&)>);
private:
	struct particle_info{
		particle_info(ParticleType type,double mass);
		ParticleType type;double mass;
		Kinematic cache;
	};
	vector<particle_info> first_vertex;
	vector<TrackProcessing> m_processing;
	double p_beam_cache;
	unsigned long m_count;
	TH1F* type_hist;
};
template <class datatype,class reaction>
class CustomAnalysis:public virtual datatype,public virtual reaction{
public:
	CustomAnalysis():Analysis(),datatype(),reaction(){}
	virtual ~CustomAnalysis(){}
protected:
	virtual void PrepairForEventAnalysis()override{
		datatype::PrepairForEventAnalysis();
	}
	virtual bool EventProcessingCondition()override{
		return datatype::EventProcessingCondition();
	}
	
	virtual bool EventPreProcessing()override{
		return reaction::EventPreProcessing();
	}
	virtual void EventPostProcessing()override{
		reaction::EventPostProcessing();
	}
	virtual bool TrackCountTrigger(int CinC,int NinC,int CinF,int NinF)override{
		return reaction::TrackCountTrigger(CinC,NinC,CinF,NinF);
	}
};
#endif
