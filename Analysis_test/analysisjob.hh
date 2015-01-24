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
enum TrackTypes {kFDN=1,kFDC=2,kCDN=11,kCDC=12};//FD neutral, FD charged, CD neutral, CD charged
enum ParticleTypes{kDummy=0,kGamma=1,kElectron=2,kPositron=3,kPi0=7,kPiPlus=8,kPiMinus=9,kNeutron=13,kProton=14,kDeuteron=45,kTriton=46,kHe3=49};
enum ForwardDetectorPlanes{kFWC1=10,kFWC2=11,kFTH1=1,kFTH2=2,kFTH3=3,kFRH1=4,kFRH2=5,kFRH3=6,kFRH4=7,kFRH5=8,kFVH=9};
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
	REventHeaderi     *fHeader;
	WTrackBank        *fTrackBankFD;
	WTrackBank        *fTrackBankCD;
	FDFTHTracks       *TrackFinderFD;
	FDEdep2Ekin       *fFDEdep2Ekin;
	WVertexBank       *fVertexBank;
	WTrackBank        *fMCTrackBank;
	WVertexBank       *fMCVertexBank;
	CCardWDET         *fDetectorTable;
	MCTrackFinder     *fMCTrackFinder;
protected:
	ClassDef(AnalysisJob,0);
};
#endif // ANALYSISJOB_H
