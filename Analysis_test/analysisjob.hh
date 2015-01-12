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
	virtual void ProcessEvent()override;
	virtual void Clear(Option_t *option = "")override;
	virtual void Print(Option_t *option = "")override;
	virtual void UserCommand(CCommand * command)override;
private:

protected:
	ClassDef(AnalysisJob,0);
};

#endif // ANALYSISJOB_H
