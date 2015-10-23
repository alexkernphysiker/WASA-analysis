// this file is distributed under 
// GPL v 3.0 license
#include <functional>
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
Analysis::Kinematic::Kinematic(){
	E=INFINITY;Th=INFINITY;Phi=INFINITY;
}
Analysis::particle_info::particle_info(ParticleType t, double m){
	type=t;mass=m;
}
void Analysis::AddParticleToFirstVertex(ParticleType type, double mass){
	first_vertex.push_back(particle_info(type,mass));
}
double Analysis::PBeam(){return p_beam_cache;}
void Analysis::CachePBeam(double value){
	if(value>0)
		p_beam_cache=value;
	else
		throw exception();
}
Analysis::Kinematic& Analysis::FromFirstVertex(ParticleType type){
	for(particle_info&info:first_vertex)
		if(info.type==type)
			return info.cache;
	throw exception();
}
void Analysis::ForFirstVertex(function<void(ParticleType,double,Kinematic&)> cyclebody){
	for(particle_info&info:first_vertex)
		cyclebody(info.type,info.mass,info.cache);
}

void Analysis::ProcessEvent(){
	m_count++;
	if(m_count%1000==0)
		Log(NoLog)<<to_string(m_count)+" events";
	SubLog log=Log(LogDebug);
	log<<"event preparing started";
	PrepairForEventAnalysis();
	log<<"event preparing finished";
	if (EventProcessingCondition()){
		log<<"event passed the condition";
		if(p_beam_cache>0){
			log<<"p_beam>0";
			if(EventPreProcessing()){
				log<<"preprocessing returned true. go further";
				int ChargedCountInCentral = fTrackBankCD->GetEntries(kCDC);
				int NeutralCountInCentral = fTrackBankCD->GetEntries(kCDN);
				int ChargedCountinForward = fTrackBankFD->GetEntries(kFDC);
				int NeutralCountinForward = fTrackBankFD->GetEntries(kFDN);
				if(TrackCountTrigger(ChargedCountInCentral,NeutralCountInCentral,ChargedCountinForward,NeutralCountinForward)){
					log<<"Track enumerating";
					vector<WTrackBank*> BANK;
					BANK.push_back(fTrackBankCD);
					BANK.push_back(fTrackBankFD);
					for(TrackProcessing&process:m_processing){
						for(WTrackBank*bank:BANK){
							WTrackIter iterator(bank);
							iterator.SetType(process.first);
							while(WTrack* track = dynamic_cast<WTrack*> (iterator.Next()))
								process.second(*track);
						}
					}
					log<<"Event postprocessing";
					EventPostProcessing();
				}log<<"Track count conditions NOT passed";
			}log<<"preprocessing returned false.";
		}else Log(LogError)<<"p_beam <= 0";
	}else log<<"event did not pass the processing condition";
}
