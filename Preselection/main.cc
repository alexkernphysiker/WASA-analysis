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
	string mode=argv[1];
	string type=argv[2];
	for(int i=2;i<=new_c;i++)
		args[i]=argv[i+2];
	gSorterConfig->ReadCmdLine(new_c,args);
	if("rec"==mode){
		if("MC_He3eta"==type)
			SetAnalysis([](){return new RE_He3eta();});
		if(
			("MC_He3pi0"==type)||
			("MC_He3pi0pi0"==type)||
			("MC_He3pi0pi0pi0"==type)
		)
			SetAnalysis([](){return new RE_He3pi0();});
	}else if("ana"==mode){
		if("MC_He3eta"==type)
			SetAnalysis([](){return new MC_He3eta();});
		if(
			("MC_He3pi0"==type)||
			("MC_He3pi0pi0"==type)||
			("MC_He3pi0pi0pi0"==type)
		)
			SetAnalysis([](){return new MC_He3pi0();});
		if("Data_He3"==type)
			SetAnalysis([](){return new Data_He3();});
	}else{
		throw math_h_error<IAnalysis>("Unknown analysis mode");
	}
	Wasa::Initialize("AnalysisModule","","RootSorter.log");
	gWasa->AddAnalysis("AnalysisModule","Raw");
	gWasa->Run();
	delete gWasa;
}
