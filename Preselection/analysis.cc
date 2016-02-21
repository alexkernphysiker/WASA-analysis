// this file is distributed under 
// MIT license
#include <functional>
#include "math_h/error.h"
#include "analysis.h"
#define static_right(A) (static_cast<decltype(A)&&>(A))
using namespace std;
IAnalysis::~IAnalysis(){}
Analysis::Analysis(){
	m_count=0;
	AddLogSubprefix("Analysis");
	TrackFinderFD = dynamic_cast<FDFTHTracks*>(gDataManager->GetAnalysisModule("FDFTHTracks","default"));
	if(TrackFinderFD!=0) fTrackBankFD = TrackFinderFD->GetTrackBank();
	CDTrackFinder = dynamic_cast<CDTracksSimple*>(gDataManager->GetAnalysisModule("CDTracksSimple","default"));
	if (CDTrackFinder!=0) fTrackBankCD = CDTrackFinder->GetTrackBank();
	fDetectorTable = dynamic_cast<CCardWDET*>(gParameterManager->GetParameterObject("CCardWDET","default")); 
	p_beam_cache=INFINITY;
	m_chain.push_back(make_pair(kFDC,TrackAnalyse::TrackProcess()));
	m_chain.push_back(make_pair(kFDN,TrackAnalyse::TrackProcess()));
	m_chain.push_back(make_pair(kCDC,TrackAnalyse::TrackProcess()));
	m_chain.push_back(make_pair(kCDN,TrackAnalyse::TrackProcess()));
}
Analysis::~Analysis(){}
//Triggers
TrackAnalyse::TrackProcess& Analysis::TrackTypeProcess(TrackType type){
	for(TrackTypeRec& rec:m_chain)
		if(type==rec.first)
			return rec.second;
	throw MathTemplates::Exception<Analysis>("Cannot find track type");
}
TrackAnalyse::EventProcess& Analysis::EventProcessing(){
	return m_eventproc;
}

bool Analysis::Trigger(int n)const{
	return DataSpecificTriggerCheck(n);
}
void Analysis::ProcessEvent(){
	m_count++;
	if(m_count%1000==0)
		Log(NoLog)<<to_string(m_count)+" events";
	SubLog log=Log(LogDebug);
	if(DataTypeSpecificEventAnalysis())
		m_eventproc.Process();
		for(const TrackTypeRec& tt:m_chain){
			vector<WTrackBank*> BANK;
			BANK.push_back(fTrackBankCD);
			BANK.push_back(fTrackBankFD);
			for(WTrackBank*bank:BANK){
				WTrackIter iterator(bank);
				while(WTrack* track = dynamic_cast<WTrack*> (iterator.Next()))
					if(track->Type()==tt.first)tt.second.Process(*track);
			}
		}
}

///BEAM MOMENTUM
double Analysis::PBeam()const{return p_beam_cache;}
void Analysis::CachePBeam(double value){
	if(value>0)p_beam_cache=value;
	else throw MathTemplates::Exception<Analysis>("Wrong beam momentum value");
}
///KINEMATICS
Analysis::Kinematic::Kinematic(){E=INFINITY;Th=INFINITY;Phi=INFINITY;}
Analysis::particle_info::particle_info(ParticleType t, double m){type=t;mass=m;}
void Analysis::AddParticleToFirstVertex(ParticleType type, double mass){
	first_vertex.push_back(particle_info(type,mass));
}
const Analysis::Kinematic& Analysis::FromFirstVertex(ParticleType type)const{
	for(const particle_info&info:first_vertex)
		if(info.type==type)
			return info.cache;
	throw MathTemplates::Exception<Analysis>("Particle not found in the vertex");
}
void Analysis::ForFirstVertex(function<void(ParticleType,double,Kinematic&)> cyclebody){
	for(particle_info&info:first_vertex)
		cyclebody(info.type,info.mass,info.cache);
}