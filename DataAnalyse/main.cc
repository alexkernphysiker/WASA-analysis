#include <Wasa.hh>
#include <SorterConfig.hh>
#include "config.h"
int main(int argc, char** argv) {
	int new_c=argc-1;
	char *args[new_c+1];
	args[0]=argv[0];
	string type=args[1];
	SetAnalysisType(type);
	for(int i=1;i<=new_c;i++)
		args[i]=argv[i+1];
	gSorterConfig->ReadCmdLine(new_c,args);
	Wasa::Initialize("AnalysisWrap","","RootSorter.log");
	gWasa->AddAnalysis("AnalysisWrap","Raw");
	gWasa->Run();
	delete gWasa;
}
