#include "gnuplot.h"
#include <math_h/functions.h>
using namespace std;
using namespace Genetic;
int main(int,char**){
	SetPlotOutput(".");
	Plot plot;
	for(double asym=-1;asym<=1;asym+=.1){
		auto F=[asym](double x){return Novosibirsk<>(x,0.0,2.0,asym);};
		plot.Function(to_string(asym),F,-2.0,2.0,0.01);
		printf("for %f:\t%f\n",asym,Sympson(F,-100.0,100.0,0.001));
	}
	return 0;
}