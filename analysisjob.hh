#ifndef ANALYSISJOB_H
#define ANALYSISJOB_H
#include <Wasa.hh>
#include <CAnalysisModule.hh>
#include "AnalyseBase/wrap.h"
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
	IAnalysis *job;
protected:
	ClassDef(AnalysisJob,0);
};
#endif // ANALYSISJOB_H
