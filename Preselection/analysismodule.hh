// this file is distributed under 
// MIT license
#ifndef MZQRVXYK
# define MZQRVXYK
#include "analysiswrap.hh"
class AnalysisModule:public AbstractAnalysis{
public:
	AnalysisModule();
	explicit AnalysisModule(const char* name);
	virtual ~AnalysisModule();
protected:
	ClassDef(AnalysisModule,1);
};

#endif