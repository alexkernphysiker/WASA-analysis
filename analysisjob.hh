#ifndef ANALYSISJOB_H
#define ANALYSISJOB_H
#include <Wasa.hh>
#include <CAnalysisModule.hh>
#include <REventWmcHeader.hh>
#include <REventHeader.hh>
#include <WTrackBank.hh>
#include <WVertexBank.hh>
#include <FDFTHTracks.hh>
#include <CDTracksSimple.hh>
#include <TH1F.h>
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
protected:
	ClassDef(AnalysisJob,0);
};
#endif // ANALYSISJOB_H
