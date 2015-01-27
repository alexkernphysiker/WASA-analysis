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
const Double_t m_target=1.875613;	//[GeV]
const Double_t m_beam=1.875613;		//[GeV]
const Double_t m_3He=2.808950;		//[GeV]
const Double_t m_p=0.938272;		//[GeV]
const Double_t m_n=0.93956;		//[GeV]
const Double_t m_4He=3.7264225;		//[GeV]
const Double_t m_eta=0.547853;		//[GeV]
enum TrackType {kFDN=1,kFDC=2,kCDN=11,kCDC=12};//FD neutral, FD charged, CD neutral, CD charged
enum ParticleType{kDummy=0,kGamma=1,kElectron=2,kPositron=3,kPi0=7,kPiPlus=8,kPiMinus=9,kNeutron=13,kProton=14,kDeuteron=45,kTriton=46,kHe3=49};
enum ForwardDetectorPlane{kFWC1=10,kFWC2=11,kFTH1=1,kFTH2=2,kFTH3=3,kFRH1=4,kFRH2=5,kFRH3=6,kFRH4=7,kFRH5=8,kFVH=9};
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
	FDEdep2Ekin	*fFDEdep2Ekin;
	WTrackBank	*fMCTrackBank;
	WVertexBank	*fMCVertexBank;
	CCardWDET	*fDetectorTable;
	MCTrackFinder	*fMCTrackFinder;
	REventWmcHeader	*fEventHeader;
protected:
	TH1F *He3_Ekin, *He3_Theta, *He3_Phi;
	ClassDef(AnalysisJob,0);
};
#endif // ANALYSISJOB_H
