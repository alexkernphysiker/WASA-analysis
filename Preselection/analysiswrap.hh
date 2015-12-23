// this file is distributed under 
// MIT license
#ifndef UWTLYUGSVUAZXUEJ
#define UWTLYUGSVUAZXUEJ
#include <Wasa.hh>
#include <CAnalysisModule.hh>
class AbstractAnalysis:CAnalysisModule{
protected:
	AbstractAnalysis();
	explicit AbstractAnalysis(const char* name,const void*data);
public:
	virtual ~AbstractAnalysis();
	virtual void ProcessEvent();
	virtual void Clear(Option_t *option = "");
	virtual void Print(Option_t *option = "");
	virtual void UserCommand(CCommand * command);
private:
	void * m_data;
};
class AnalysisModule:public AbstractAnalysis{
public:
	AnalysisModule();
	explicit AnalysisModule(const char* name);
	virtual ~AnalysisModule();
protected:
	ClassDef(AnalysisModule,0);
};
class ReconstructionModule:public AbstractAnalysis{
public:
	ReconstructionModule();
	explicit ReconstructionModule(const char* name);
	virtual ~ReconstructionModule();
protected:
	ClassDef(ReconstructionModule,1);
};
#endif
