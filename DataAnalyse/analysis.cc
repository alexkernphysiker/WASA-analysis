#include <functional>
#include "analysis.h"
using namespace std;
IAnalysis::~IAnalysis(){}
Analysis::Analysis():Logger("Analysis"){
	LogMessage(LogDebug,"Constructor called");
	checkprepared=false;
	TrackFinderFD = dynamic_cast<FDFTHTracks*>(gDataManager->GetAnalysisModule("FDFTHTracks","default"));
	if(TrackFinderFD!=0) fTrackBankFD = TrackFinderFD->GetTrackBank();
	CDTrackFinder = dynamic_cast<CDTracksSimple*>(gDataManager->GetAnalysisModule("CDTracksSimple","default"));
	if (CDTrackFinder!=0) fTrackBankCD = CDTrackFinder->GetTrackBank();
	fDetectorTable = dynamic_cast<CCardWDET*>(gParameterManager->GetParameterObject("CCardWDET","default")); 
}
Analysis::~Analysis(){
	LogMessage(LogDebug,"Destructor called");
}
typedef pair<WTrackBank*,function<bool(WTrack*,TVector3)>> DetectorToProcess;
void Analysis::ProcessEvent(){
	LogMessage(LogDebug,"==== begin processing event ====");
	if(!checkprepared){
		PrepareCheck();
		checkprepared=true;
	}
	if (EventProcessingCondition()){
		LogMessage(LogDebug,"event passed the condition");
		double beam_momenta=PBeam();
		if(beam_momenta>0){
			LogMessage(LogDebug,"p_beam>0");
			TVector3 p_beam;
			p_beam.SetMagThetaPhi(beam_momenta,0,0);
			if(EventPreProcessing(p_beam)){
				LogMessage(LogDebug,"preprocessing returned true. go further");
				int ChargedCountInCentral = fTrackBankCD->GetEntries(kCDC);
				int NeutralCountInCentral = fTrackBankCD->GetEntries(kCDN);
				int ChargedCountinForward = fTrackBankFD->GetEntries(kFDC);
				if(TrackCountTrigger(ChargedCountInCentral,NeutralCountInCentral,ChargedCountinForward)){
					LogMessage(LogDebug,"Track count conditions passed");
					DetectorToProcess CENTRAL=make_pair(fTrackBankCD,[this](WTrack* t,TVector3 p){
						LogMessage(LogDebug,"CENTRAL tracks processing");
						return CentralTrackProcessing(t,p);
					});
					DetectorToProcess FORWARD=make_pair(fTrackBankFD,[this](WTrack* t,TVector3 p){
						LogMessage(LogDebug,"FORWARD tracks processing");
						return ForwardTrackProcessing(t,p);
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
						for (int trackindex=0; trackindex<NrTracks; trackindex++) { 
							auto track=DETECTOR.first->GetTrack(trackindex);
							if(!DETECTOR.second(track,p_beam))
								return;
						}
					}
					LogMessage("Event postprocessing");
					EventPostProcessing(p_beam);
				}LogMessage(LogDebug,"Track count conditions NOT passed");
			}LogMessage(LogDebug,"preprocessing returned false.");
		}LogMessage(LogWarning,"p_beam <= 0");
	}else LogMessage(LogDebug,"event did not pass the condition");
	LogMessage(LogDebug,"==== end processing event ====");
}
