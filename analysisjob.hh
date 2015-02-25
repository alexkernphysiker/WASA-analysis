#ifndef ANALYSISJOB_H
#define ANALYSISJOB_H
#include "CAnalysisModule.hh"
#include "REventWmcHeader.hh"
#include "REventHeader.hh"
#include "FDFTHTracks.hh"
#include "CDTracksSimple.hh"
#include "FDEdep2Ekin.hh"
#include "CCardWDET.hh"
#include "WTrackBank.hh"
#include "WVertexBank.hh"
#include <fstream>
#include <TH1.h>
#include <TH2.h>
#include "TCutG.h"
class AnalysisJob : public CAnalysisModule{
public:
	AnalysisJob();
	explicit AnalysisJob(const char * name);
	virtual ~AnalysisJob();
	virtual void ProcessEvent();
	virtual void Clear(Option_t *option = "");
	virtual void Print(Option_t *option = "");
	virtual void UserCommand(CCommand * command);
private:
	REventHeader	*fHeader;
	WTrackBank	*fTrackBankFD;
	WTrackBank	*fTrackBankCD;
	FDFTHTracks	*TrackFinderFD;
	CDTracksSimple	*TrackFinderCD;
	WTrackBank	*fMCTrackBank;
	WVertexBank	*fMCVertexBank;
	CCardWDET	*fDetectorTable;
	MCTrackFinder	*fMCTrackFinder;
	REventWmcHeader	*fEventHeader;
protected:
	ClassDef(AnalysisJob,0);
};
#endif // ANALYSISJOB_H
