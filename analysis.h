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
protected:
	Analysis();
public:
	virtual ~Analysis();
	virtual void ProcessEvent()override;
protected:
	virtual bool ProcessingCondition()=0;
	virtual double PBeam()=0;
	virtual void SpecificProcessing()=0;
	virtual double EventWeight()=0;
	WTrackBank *fTrackBankFD,*fTrackBankCD;
	//FDEdep2Ekin *He3DepKin;
	CCardWDET *fDetectorTable;
	TH1F *P_Beam;
};
class MCAnalysis:public virtual Analysis{
public:
    MCAnalysis();
    virtual ~MCAnalysis();
protected:
	virtual bool ProcessingCondition()override;
	virtual double PBeam()override;
	virtual void SpecificProcessing()override;
	virtual double EventWeight()override;
	WTrackBank *fMCTrackBank;
	WVertexBank *fMCVertexBank;
	REventWmcHeader   *fEventHeader;
	TH1F *He3_Ekin,*He3_Theta,*He3_Phi;
};
const double m_p=0.938272;//[GeV]
const double m_n=0.93956;//[GeV]
const double m_d=1.875613;//[GeV]
const double m_3He=2.808950;//[GeV]
const double m_4He=3.7264225;//[GeV]
const double m_eta=0.547853;//[GeV]
enum TrackType{kFDN=1,kFDC=2,kCDN=11,kCDC=12};
enum ForwardDetectorPlane{kFWC1=10,kFWC2=11,kFTH1=1,kFTH2=2,kFTH3=3,kFRH1=4,kFRH2=5,kFRH3=6,kFRH4=7,kFRH5=8,kFVH=9};

enum ParticleType{
	kDummy=0,kGamma=1,kElectron=2,kPositron=3,kPi0=7,kPiPlus=8,kPiMinus=9,
	kNeutron=13,kProton=14,kEta=17,kDeuteron=45,kTriton=46,kHe3=49
};
#endif
