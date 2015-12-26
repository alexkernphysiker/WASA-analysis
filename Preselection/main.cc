// this file is distributed under 
// MIT license
#include <Wasa.hh>
#include <SorterConfig.hh>
#include "config.h"
#include "he3.h"
using namespace std;
int main(int argc, char** argv) {
	InitLog(LogWarning,argv[1]);
	int new_c=argc-2;
	char *args[new_c+1];
	args[0]=argv[0];
	SetAnalysisType(argv[1],argv[2]);
	for(int i=2;i<=new_c;i++)
		args[i]=argv[i+2];
	gSorterConfig->ReadCmdLine(new_c,args);
	Wasa::Initialize("AnalysisModule","","RootSorter.log");
	gWasa->AddAnalysis("AnalysisModule","Raw");
	gWasa->Run();
	delete gWasa;
}
