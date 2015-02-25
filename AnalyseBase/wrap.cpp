#include <memory>
#include "wrap.h"
#include "analysis.h"
using namespace std;
IAnalysis::~IAnalysis(){}
IAnalysis* CreateAnalysis(){
	return new Analysis();
}
