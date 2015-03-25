#ifndef ____PRCbZGGy
#define ____PRCbZGGy
#include <Wasa.hh>
#include <CAnalysisModule.hh>
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
	void * m_data;
protected:
	ClassDef(AnalysisJob,0);
};
#endif
