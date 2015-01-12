#include "Wasa.hh"
#include "SorterConfig.hh"

int main(int argc, char** argv) {
  gSorterConfig->ReadCmdLine(argc,argv);
  Wasa::Initialize("examplejob","","RootSorter.log");
  gWasa->AddAnalysis("examplejob","Raw");
  gWasa->Run();
  delete gWasa;
}
