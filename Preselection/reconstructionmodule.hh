// this file is distributed under 
// MIT license
#ifndef QYSGRPEP
# define QYSGRPEP
#include "abstractanalysis.hh"
class ReconstructionModule:public AbstractAnalysis{
public:
	ReconstructionModule();
	explicit ReconstructionModule(const char* name);
	virtual ~ReconstructionModule();
protected:
	ClassDef(ReconstructionModule,2);
};
#endif