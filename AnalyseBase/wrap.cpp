#include <memory>
#include "wrap.h"
#include "analysis.h"
using namespace std;
IAnalysis* CreateAnalysis(){
	return new Analysis();
}
