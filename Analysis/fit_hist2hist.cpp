// this file is distributed under 
// GPL v 3.0 license
#include "fit_hist2hist.h"
using namespace std;
using namespace Genetic;
Hist2Hist::Hist2Hist(shared_ptr< hist > data, shared_ptr< hist > mc, Hist2Hist::bg_func background){
	DATA=data;MC=mc;BG=background;
}
Hist2Hist::~Hist2Hist(){}
hist::point&Get(std::shared_ptr<hist>H,size_t i,double corr){
	int newi=int(corr+i);
	if(newi<0)return H->operator[](0);
	size_t last=H->count()-1;
	if(newi>last)return H->operator[](last);
	return H->operator[](newi);
}
double Hist2Hist::operator()(ParamSet&& P){
	if(P[0]<0)throw exception();
	if(DATA->count()!=MC->count())
		throw exception();
	double res=0;
	double z=DATA->count()-P.Count();
	if(z<=0)throw exception();
	for(size_t i=0,n=DATA->count();i<n;i++){
		double x=MC->operator[](i).x;
		if(x!=DATA->operator[](i).x)throw exception();
		double y=Get(MC,i,P[1]).y+BG(x,static_right(P)),dy=Get(MC,i,P[1]).dy,
			Y=DATA->operator[](i).y,dY=DATA->operator[](i).dy;
		y*=P[0];dy*=P[0];
		if((dy<=0)||(dY<=0))throw exception();
		res+=pow((Y-y)/(dY+dy),2);
	}
	return res/z;
}
