// this file is distributed under 
// GPL v 3.0 license
#include "fit_hist2hist.h"
using namespace std;
using namespace Genetic;
Hist2Hist::Hist2Hist(shared_ptr< hist > data, shared_ptr< hist > mc, Hist2Hist::bg_func background){
	DATA=data;MC=mc;BG=background;
}
Hist2Hist::~Hist2Hist(){}
double Hist2Hist::operator()(ParamSet&& P){
	if(DATA->count()!=MC->count())
		throw exception();
	double res=0;
	double z=DATA->count()-P.Count();
	if(z<=0)throw exception();
	for(size_t i=0,n=DATA->count();i<n;i++){
		double x=MC->operator[](i).x;
		if(x!=DATA->operator[](i).x)throw exception();
		double y=MC->operator[](i).y+BG(x,static_right(P)),dy=MC->operator[](i).dy,
			Y=DATA->operator[](i).y,dY=DATA->operator[](i).dy;
		y*=P[0];dy*=P[0];
		res+=pow((Y-y)/(dY+dy),2);
	}
	return res/z;
}
