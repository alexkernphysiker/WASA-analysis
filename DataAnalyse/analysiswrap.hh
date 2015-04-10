#ifndef UWTLYUGSVUAZXUEJ
#define UWTLYUGSVUAZXUEJ
#include <Wasa.hh>
#include <CAnalysisModule.hh>
void SetAnalysisType(char *type);
class AnalysisWrap:public CAnalysisModule{
public:
	AnalysisWrap();
	explicit AnalysisWrap(const char* name);
	virtual ~AnalysisWrap();
	virtual void ProcessEvent();
	virtual void Clear(Option_t *option = "");
	virtual void Print(Option_t *option = "");
	virtual void UserCommand(CCommand * command);
private:
	void * m_data;
protected:
	ClassDef(AnalysisWrap,0);
};
#endif
