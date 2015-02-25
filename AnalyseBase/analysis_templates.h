#ifndef ___ANALYSIS_TEMPLATES_______
#define ___ANALYSIS_TEMPLATES_______
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
template<>
class Analysis : public CAnalysisModule{
public:
	Analysis();
	explicit Analysis(const char * name);
	virtual ~Analysis();
	virtual void ProcessEvent();
	virtual void Clear(Option_t *option = "");
	virtual void Print(Option_t *option = "");
	virtual void UserCommand(CCommand * command);
private:
	REventHeader	*fHeader;
	WTrackBank	*fTrackBankFD;
	WTrackBank	*fTrackBankCD;
};
#endif