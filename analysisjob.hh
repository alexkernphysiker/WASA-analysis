#ifndef ____PRCbZGGy
#define ____PRCbZGGy
#include <Wasa.hh>
#include <CAnalysisModule.hh>
class MCHe3Eta : public CAnalysisModule{
public:
	MCHe3Eta();
	explicit MCHe3Eta(const char * name);
	virtual ~MCHe3Eta();
	virtual void ProcessEvent();
	virtual void Clear(Option_t *option = "");
	virtual void Print(Option_t *option = "");
	virtual void UserCommand(CCommand * command);
private:
	void * m_data;
protected:
	ClassDef(MCHe3Eta,0);
};
#endif
