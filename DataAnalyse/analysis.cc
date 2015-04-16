#include <functional>
#include "analysis.h"
using namespace std;
IAnalysis::~IAnalysis(){}
Analysis::Analysis(){
	checkprepared=false;
	TrackFinderFD = dynamic_cast<FDFTHTracks*>(gDataManager->GetAnalysisModule("FDFTHTracks","default"));
	if(TrackFinderFD!=0) fTrackBankFD = TrackFinderFD->GetTrackBank();
	CDTrackFinder = dynamic_cast<CDTracksSimple*>(gDataManager->GetAnalysisModule("CDTracksSimple","default"));
	if (CDTrackFinder!=0) fTrackBankCD = CDTrackFinder->GetTrackBank();
	fDetectorTable = dynamic_cast<CCardWDET*>(gParameterManager->GetParameterObject("CCardWDET","default")); 
}
Analysis::~Analysis(){}
typedef pair<WTrackBank*,function<bool(WTrack*,TVector3)>> DetectorToProcess;
void Analysis::ProcessEvent(){
	if(!checkprepared){
		PrepareCheck();
		checkprepared=true;
	}
	if (EventProcessingCondition()){
		double beam_momenta=PBeam();
		if(beam_momenta>0){
			TVector3 p_beam;
			p_beam.SetMagThetaPhi(beam_momenta,0,0);
			if(EventPreProcessing(p_beam)){
				int ChargedCountInCentral = fTrackBankCD->GetEntries(kCDC);
				int NeutralCountInCentral = fTrackBankCD->GetEntries(kCDN);
				int ChargedCountinForward = fTrackBankFD->GetEntries(kFDC);
				if(TrackCountTrigger(ChargedCountInCentral,NeutralCountInCentral,ChargedCountinForward)){
					DetectorToProcess CENTRAL=make_pair(fTrackBankCD,[this](WTrack* t,TVector3 p){return CentralTrackProcessing(t,p);});
					DetectorToProcess FORWARD=make_pair(fTrackBankFD,[this](WTrack* t,TVector3 p){return ForwardTrackProcessing(t,p);});
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
					EventPostProcessing(p_beam);
				}
			}
		}
	}
}
