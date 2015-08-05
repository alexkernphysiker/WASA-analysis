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
}
Analysis::~Analysis(){}
typedef pair<WTrackBank*,function<bool(WTrack&&,TVector3)>> DetectorToProcess;
void Analysis::ProcessEvent(){
	m_count++;
	if(m_count%1000==0)
		Log(NoLog)<<to_string(m_count)+" events";
	SubLog log=Log(LogDebug);
	if (EventProcessingCondition()){
		log<<"event passed the condition";
		double beam_momenta=PBeam();
		if(beam_momenta>0){
			log<<"p_beam>0";
			TVector3 p_beam;
			p_beam.SetMagThetaPhi(beam_momenta,0,0);
			if(EventPreProcessing(static_right(p_beam))){
				log<<"preprocessing returned true. go further";
				int ChargedCountInCentral = fTrackBankCD->GetEntries(kCDC);
				int NeutralCountInCentral = fTrackBankCD->GetEntries(kCDN);
				int ChargedCountinForward = fTrackBankFD->GetEntries(kFDC);
				if(TrackCountTrigger(ChargedCountInCentral,NeutralCountInCentral,ChargedCountinForward)){
					log<<"Track count conditions passed";
					DetectorToProcess CENTRAL=make_pair(fTrackBankCD,[this,&log](WTrack&& t,TVector3&& p){
						log<<"CENTRAL tracks processing";
						return CentralTrackProcessing(static_right(t),static_right(p));
					});
					DetectorToProcess FORWARD=make_pair(fTrackBankFD,[this,&log](WTrack&& t,TVector3&& p){
						log<<"FORWARD tracks processing";
						return ForwardTrackProcessing(static_right(t),static_right(p));
					});
					vector<DetectorToProcess> QUEUE;
					if(CentralFirst()){
						QUEUE.push_back(CENTRAL);
						QUEUE.push_back(FORWARD);
					}else{
						QUEUE.push_back(FORWARD);
						QUEUE.push_back(CENTRAL);
					}
					for(DetectorToProcess DETECTOR:QUEUE){
						int NrTracks=DETECTOR.first->GetEntries();
						for (int trackindex=0; trackindex<NrTracks; trackindex++)
							if(!DETECTOR.second(static_cast<WTrack&&>(*DETECTOR.first->GetTrack(trackindex)),static_right(p_beam))){
								log<<"proceccing returned FALSE. Exiting event processing.";
								return;
							}
					}
					log<<"Event postprocessing";
					EventPostProcessing(static_right(p_beam));
				}log<<"Track count conditions NOT passed";
			}log<<"preprocessing returned false.";
		}else Log(LogError)<<"p_beam <= 0";
	}else log<<"event did not pass the condition";
}
