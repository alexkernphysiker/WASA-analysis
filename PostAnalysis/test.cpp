#include "gnuplot.h"
#include <math_h/functions.h>
using namespace std;
using namespace Genetic;
int main(int,char**){
	SetPlotOutput(".");
	Plot plot;
	for(double asym=-1;asym<=1;asym+=0.1){
		auto F=[asym](double x){return Novosibirsk<>(x,0.0,1.0,asym);};
		plot.Function(to_string(asym),F,-10,10,0.01);
		printf("for %f:\t%f\n",asym,Sympson(F,-50.0,50.0,0.01));
	}
	return 0;
}