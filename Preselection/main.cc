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
		//ToDo: provide reconstruction algorithms
	}else if("ana"==mode){
		IAnalysis *alg=nullptr;
		if("MC_He3eta"==type)alg=new MC_He3eta();
		if(("MC_He3pi0"==type)||("MC_He3pi0pi0"==type)||("MC_He3pi0pi0pi0"==type))alg=new MC_He3pi0();
		if("Data_He3"==type)alg=new Data_He3();
		SetAnalysis(alg);
	}else{
		throw math_h_error<IAnalysis>("Unknown analysis mode");
	}
	Wasa::Initialize("AnalysisModule","","RootSorter.log");
	gWasa->AddAnalysis("AnalysisModule","Raw");
	gWasa->Run();
	delete gWasa;
}
