#ifndef UWTLYUGSVUAZXUEJ
#define UWTLYUGSVUAZXUEJ
#include <Wasa.hh>
#include <CAnalysisModule.hh>
class AnalysisWrap:public CAnalysisModule{
protected:
    AnalysisWrap();
    explicit AnalysisWrap(const char* name);
public:
	virtual ~AnalysisWrap();
	virtual void ProcessEvent();
	virtual void Clear(Option_t *option = "");
	virtual void Print(Option_t *option = "");
	virtual void UserCommand(CCommand * command);
protected:
	void * m_data;
};
//particular types of analysis 
class MCHe3Eta : public AnalysisWrap{
public:
	MCHe3Eta();
	explicit MCHe3Eta(const char * name);
	virtual ~MCHe3Eta();
protected:
	ClassDef(MCHe3Eta,0);
};
class DataHe3Eta : public AnalysisWrap{
public:
	DataHe3Eta();
	explicit DataHe3Eta(const char * name);
	virtual ~DataHe3Eta();
protected:
	ClassDef(DataHe3Eta,1);
};
#endif
