#ifndef ROUTINES_H
#define ROUTINES_H
#include <Wasa.hh>
#include <CDataManager.hh>
#include <CHistoManager.hh>
class IAnalysis{
public:
	virtual ~IAnalysis();
	virtual Init(CDataManager *dataManager, CHistoManager *histoManager)=0;
	virtual Processevent(Wasa* wasa)=0;
};
IAnalysis* CreateAnalysis();
#endif // ROUTINES_H
