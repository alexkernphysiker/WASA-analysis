#include <Wasa.hh>
#include <SorterConfig.hh>
#include "config.h"
using namespace std;
int main(int argc, char** argv) {
	InitLog(LogWarning,argv[1]);
	int new_c=argc-1;
	char *args[new_c+1];
	args[0]=argv[0];
	SetAnalysisType(argv[1]);
	for(int i=1;i<=new_c;i++)
		args[i]=argv[i+1];
	gSorterConfig->ReadCmdLine(new_c,args);
	Wasa::Initialize("AnalysisWrap","","RootSorter.log");
	gWasa->AddAnalysis("AnalysisWrap","Raw");
	gWasa->Run();
	delete gWasa;
}
