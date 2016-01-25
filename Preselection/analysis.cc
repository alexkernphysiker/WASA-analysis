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
}
Analysis::~Analysis(){}
//Triggers
Analysis::TriggerProcess::TriggerProcess(int n):N(n){}
Analysis::TriggerProcess::~TriggerProcess(){}
Analysis::TriggerProcess::const_iterator Analysis::TriggerProcess::begin()const{
	return m_chain.begin();
}
Analysis::TriggerProcess::const_iterator Analysis::TriggerProcess::cbegin() const{
	return m_chain.cbegin();
}
Analysis::TriggerProcess::const_iterator Analysis::TriggerProcess::end() const{
	return m_chain.end();
}
Analysis::TriggerProcess::const_iterator Analysis::TriggerProcess::cend() const{
	return m_chain.cend();
}
size_t Analysis::TriggerProcess::size() const{
	return m_chain.size();
}
int Analysis::TriggerProcess::number() const{
	return N;
}
TrackAnalyse::TrackProcess& Analysis::TriggerProcess::TrackTypeProcess(TrackType type){
	for(TrackTypeRec& rec:m_chain)
		if(type==rec.first)
			return rec.second;
	m_chain.push_back(make_pair(type,TrackAnalyse::TrackProcess()));
	for(TrackTypeRec& rec:m_chain)
		if(type==rec.first)
			return rec.second;
	throw MathTemplates::Exception<TriggerProcess>("Cannot add track type");
}
Analysis::TriggerProcess& Analysis::Trigger(int n){
	for(TriggerProcess& trigger:m_triggers)
		if(trigger.number()==n)
			return trigger;
	m_triggers.push_back(TriggerProcess(n));
	for(TriggerProcess& trigger:m_triggers)
		if(trigger.number()==n)
			return trigger;
	throw MathTemplates::Exception<Analysis>("Cannot add trigger");
}

void Analysis::ProcessEvent(){
	m_count++;
	if(m_count%1000==0)
		Log(NoLog)<<to_string(m_count)+" events";
	SubLog log=Log(LogDebug);
	if(DataTypeSpecificEventAnalysis()){
		for(const TriggerProcess& trigger:m_triggers)
			if(DataSpecificTriggerCheck(trigger.number())){
				vector<WTrackBank*> BANK;
				BANK.push_back(fTrackBankCD);
				BANK.push_back(fTrackBankFD);
				for(WTrackBank*bank:BANK){
					WTrackIter iterator(bank);
					while(WTrack* track = dynamic_cast<WTrack*> (iterator.Next())){
						for(const TriggerProcess::TrackTypeRec&TT:trigger){
							if(track->Type()==TT.first)
								TT.second.Process(*track);
						}
					}
				}
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
Analysis::Kinematic& Analysis::FromFirstVertex(ParticleType type){
	for(particle_info&info:first_vertex)
		if(info.type==type)
			return info.cache;
	throw MathTemplates::Exception<Analysis>("Particle not found in the vertex");
}
void Analysis::ForFirstVertex(function<void(ParticleType,double,Kinematic&)> cyclebody){
	for(particle_info&info:first_vertex)
		cyclebody(info.type,info.mass,info.cache);
}