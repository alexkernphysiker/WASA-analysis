#ifndef _____EMyyUotT
#define _____EMyyUotT
#include <memory>
#include <TLorentzVector.h>
#include <CAnalysisModule.hh>
#include <CDataManager.hh>
#include <FDEdep2Ekin.hh>
#include <CCardWDET.hh>
#include "math_h/sigma.h"
#include "analysisjob.hh"
class IAnalysis{
public:
	virtual ~IAnalysis();
	virtual void ProcessEvent()=0;
};
class Analysis:public virtual IAnalysis{
public:
	Analysis();
	virtual ~Analysis();
	virtual void ProcessEvent()override;
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
const Double_t m_p=0.938272;//[GeV]
const Double_t m_n=0.93956;//[GeV]
const Double_t m_d=1.875613;//[GeV]
const Double_t m_3He=2.808950;//[GeV]
const Double_t m_4He=3.7264225;//[GeV]
const Double_t m_eta=0.547853;//[GeV]
enum TrackType{kFDN=1,kFDC=2,kCDN=11,kCDC=12};
enum ParticleType{kDummy=0,kGamma=1,kElectron=2,kPositron=3,kPi0=7,kPiPlus=8,kPiMinus=9,kNeutron=13,kProton=14,kDeuteron=45,kTriton=46,kHe3=49};
enum ForwardDetectorPlane{kFWC1=10,kFWC2=11,kFTH1=1,kFTH2=2,kFTH3=3,kFRH1=4,kFRH2=5,kFRH3=6,kFRH4=7,kFRH5=8,kFVH=9};
#endif