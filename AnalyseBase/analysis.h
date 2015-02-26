#ifndef ANALYSIS_H
#define ANALYSIS_H
#include "wrap.h"
#include <REventWmcHeader.hh>
#include <REventHeader.hh>
#include <WTrackBank.hh>
#include <WVertexBank.hh>
#include <FDFTHTracks.hh>
#include <CDTracksSimple.hh>
#include <TH1F.h>
class Analysis : public IAnalysis{
public:
	Analysis();
	virtual ~Analysis();
	virtual void Init(CDataManager *dataManager, CHistoManager *histoManager)override;
	virtual void Processevent(Wasa* wasa) override;
private:
	REventHeader	*fHeader;
	REventWmcHeader	*fEventHeader;
	WTrackBank	*fTrackBankFD;
	WTrackBank	*fTrackBankCD;
	FDFTHTracks	*TrackFinderFD;
	CDTracksSimple	*TrackFinderCD;
	WTrackBank	*fMCTrackBank;
	WVertexBank	*fMCVertexBank;
	CCardWDET	*fDetectorTable;
	MCTrackFinder	*fMCTrackFinder;
	TH1F *He3_Ekin, *He3_Theta, *He3_Phi;
};

#endif // ANALYSIS_H
