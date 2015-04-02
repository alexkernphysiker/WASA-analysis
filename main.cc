#include "Wasa.hh"
#include "SorterConfig.hh"
int main(int argc, char** argv) {
	int new_c=argc-2;
	char *args[new_c+1];
	args[0]=argv[0];
	for(int i=1;i<=new_c;i++)
		args[i]=argv[i+2];
	gSorterConfig->ReadCmdLine(new_c,args);
	Wasa::Initialize(argv[1],"","RootSorter.log");
	gWasa->AddAnalysis(argv[1],argv[2]);
	gWasa->Run();
	delete gWasa;
}
