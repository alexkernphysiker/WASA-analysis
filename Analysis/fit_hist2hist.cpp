// this file is distributed under 
// GPL v 3.0 license
#include "fit_hist2hist.h"
using namespace std;
using namespace Genetic;

hist::point&Get(std::shared_ptr<hist>H,size_t i,double corr){
	int newi=int(-corr+i);
	if(newi<0)return H->operator[](0);
	size_t last=H->count()-1;
	if(newi>last)return H->operator[](last);
	return H->operator[](newi);
}
HistToHists::HistToHists(shared_ptr<hist> data, vector<shared_ptr<hist>> simulation):m_simulation(simulation){
	m_data=data;
}
HistToHists::~HistToHists(){}
double HistToHists::operator()(const ParamSet& P)const{
	if(P.Count()!=(m_simulation.size()+1))
		throw exception();
	for(auto MC:m_simulation)
		if(MC->count()!=m_data->count())throw exception();
	double res=0;
	double z=m_data->count()-P.Count();
	if(z<=0)throw exception();
	for(size_t i=0,n=m_data->count();i<n;i++){
		double x=m_data->operator[](i).x;
		double y,dy;
		size_t index=1;
		for(auto MC:m_simulation){
			if(x!=MC->operator[](i).x)throw exception();
			y=Get(MC,i,P[0]).y*P[index];
			dy=Get(MC,i,P[0]).dy*P[index];
			index++;
		}
		double Y=m_data->operator[](i).y,dY=m_data->operator[](i).dy;
		y*=P[1];dy*=P[1];
		if((dy<0)||(dY<=0))throw exception();
		res+=pow((Y-y)/(dY+dy),2);
	}
	return res/z;
}
shared_ptr<hist> HistToHists::GetSubHist(size_t i, ParamSet&& P){
	if(P.Count()!=(m_simulation.size()+1))
		throw exception();
	if(P[i+1]<0)
		throw exception();
	if(i>=m_simulation.size())
		throw exception();
	auto res=make_shared<hist>(*(m_simulation[i]));
	res->operator*=(P[i+1]);
	int ofs=int(P[0]);
	if(ofs<0)res->operator<<(-ofs);
	if(ofs>0)res->operator>>(ofs);
	return res;
}

