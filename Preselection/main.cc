// this file is distributed under 
// MIT license
#include <Wasa.hh>
#include <SorterConfig.hh>
#include "config.h"
using namespace std;
int main(int argc, char** argv) {
	InitLog(LogWarning,argv[1]);
	int new_c=argc-2;
	char *args[new_c+1];
	args[0]=argv[0];
	string mode=argv[1];
	for(int i=2;i<=new_c;i++)
		args[i]=argv[i+2];
	gSorterConfig->ReadCmdLine(new_c,args);
	if("rec"==mode){
		SetReconstructionType(argv[2]);
		Wasa::Initialize("ReconstructionModule","","RootSorter.log");
		gWasa->AddAnalysis("ReconstructionModule","Raw");
	}
	if("ana"==mode){
		SetAnalysisType(argv[2]);
		Wasa::Initialize("AnalysisModule","","RootSorter.log");
		gWasa->AddAnalysis("AnalysisModule","Raw");
	}
	gWasa->Run();
	delete gWasa;
}
