// this file is distributed under 
// MIT license
#include <fstream>
#include <math_h/exception_math_h.h>
#include <str_get.h>
#include "reconstruction_fit.h"
using namespace std;
using namespace Genetic;
string SimulationDataPath(){
	static string str="";
	if(str=="")str=ENV(PRESEL_DATA)+string("/../Reconstruction/");
	return str;
}
shared_ptr<FitPoints> WeightPoints(
	shared_ptr<FitPoints>src,
	ParamSet&&pos,ParamSet&&delta,
	size_t index
){
	auto result=make_shared<FitPoints>();
	if(pos.size()!=index)
		throw math_h_error<FitPoints>("weighting error");
	if(delta.size()!=index)
		throw math_h_error<FitPoints>("weighting error");
	if(index>src->dimensions())
		throw math_h_error<FitPoints>("weighting error");
	if(index==src->dimensions()){
		double beg=src->Ymin(),end=src->Ymax(),d=(end-beg)/50.0;
		for(double y=beg;y<=end;y+=d){
			double weight=0;
			for(Point&P:(*src)){
				bool status=true;
				for(size_t i=0;(i<delta.size())&&status;i++)status=(pow(pos[i]-P.X()[i],2)<pow(delta[i],2));
				if(status)weight+=1.0;
			}
			result<<Point(pos,y,weight);
		}
	}else{
		double beg=src->min()[index],end=src->max()[index],d=(end-beg)/50.0;
		ParamSet new_pos=pos,new_delta=delta;
		new_pos<<0;delta<<d;
		for(double x=beg;x<=end;x+=d){
			new_pos[index]=x;
			result<<WeightPoints(src,static_cast<ParamSet&&>(pos),static_cast<ParamSet&&>(delta),index+1);
		}
	}
	return result;
}

