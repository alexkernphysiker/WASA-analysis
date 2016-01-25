// this file is distributed under 
// MIT license
#ifndef QTBEVPHGEADUWGMY
#define QTBEVPHGEADUWGMY
#include <memory>
#include <utility>
#include <vector>
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
#include "trackprocessing.h"
#include "analysismodule.hh"
#include "log.h"
enum TrackType{kFDN=1,kFDC=2,kCDN=11,kCDC=12};
enum ParticleType{
	kDummy=0,kGamma=1,kElectron=2,kPositron=3,kPi0=7,kPiPlus=8,kPiMinus=9,
	kNeutron=13,kProton=14,kEta=17,kDeuteron=45,kTriton=46,kHe3=49
};
class Analysis:public virtual IAnalysis,public virtual Logger{
protected:
	Analysis();
public:
	virtual ~Analysis();
	virtual void ProcessEvent()final;
protected:
	class TriggerProcess{
	public:
		TriggerProcess(int n);
		~TriggerProcess();
		int number()const;
		size_t size()const;
		TrackAnalyse::TrackProcess&TrackTypeProcess(TrackType);
		typedef std::pair<TrackType,TrackAnalyse::TrackProcess> TrackTypeRec;
		typedef std::vector<TrackTypeRec> TrackTypeRecs;
		typedef TrackTypeRecs::const_iterator const_iterator;
		const_iterator begin()const;
		const_iterator cbegin()const;
		const_iterator end() const;
		const_iterator cend() const;
	private:
		int N;
		TrackTypeRecs m_chain;
	};
	TriggerProcess&Trigger(int n);
	virtual bool DataTypeSpecificEventAnalysis()=0;
	virtual bool DataSpecificTriggerCheck(int n)=0;
	//Beam momenta calculation
	double PBeam()const;
	void CachePBeam(double value);
	//Kinematics calculations
	void AddParticleToFirstVertex(ParticleType type,double mass);
	struct Kinematic{Kinematic();double E,Th,Phi;};
	Kinematic&FromFirstVertex(ParticleType type);
	void ForFirstVertex(std::function<void(ParticleType,double,Kinematic&)>);
private:
	FDFTHTracks* TrackFinderFD;
	CDTracksSimple* CDTrackFinder;
	WTrackBank *fTrackBankFD,*fTrackBankCD;
	CCardWDET *fDetectorTable;
	std::vector<TriggerProcess> m_triggers;
	struct particle_info{
		particle_info(ParticleType type,double mass);
		ParticleType type;double mass;
		Kinematic cache;
	};
	vector<particle_info> first_vertex;
	double p_beam_cache;
	unsigned long m_count;
};
template <class datatype,class reaction>
class CustomAnalysis:public virtual datatype,public virtual reaction{
public:
	CustomAnalysis():Analysis(),datatype(),reaction(){}
	virtual ~CustomAnalysis(){}
protected:
	virtual bool DataTypeSpecificEventAnalysis()override{
		return datatype::DataTypeSpecificEventAnalysis();
	}
	virtual bool DataSpecificTriggerCheck(int n)override{
		return datatype::DataSpecificTriggerCheck(n);
	}
};
#endif
