#include "Wasa.hh"
#include "SorterConfig.hh"

int main(int argc, char** argv) {
  gSorterConfig->ReadCmdLine(argc,argv);
  Wasa::Initialize("AnalysisJob","","RootSorter.log");
  gWasa->AddAnalysis("AnalysisJob","Raw");
  gWasa->Run();
  delete gWasa;
}
