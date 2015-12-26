// this file is distributed under 
// MIT license
#ifndef MZQRVXYK
# define MZQRVXYK
#include <Wasa.hh>
#include <CAnalysisModule.hh>
class IAnalysis{
public:
	virtual ~IAnalysis();
	virtual void ProcessEvent()=0;
};
class AnalysisModule:public CAnalysisModule{
public:
	AnalysisModule();
	explicit AnalysisModule(const char* name);
	virtual ~AnalysisModule();
	virtual void ProcessEvent();
	virtual void Clear(Option_t *option = "");
	virtual void Print(Option_t *option = "");
	virtual void UserCommand(CCommand * command);
private:
	IAnalysis*m_data;
protected:
	ClassDef(AnalysisModule,0);
};

#endif