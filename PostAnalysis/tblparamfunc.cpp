#include "tblparamfunc.h"
using namespace Genetic;
TblParamFunc::TblParamFunc(FuncTbl func){
	m_function=func;
}
TblParamFunc::~TblParamFunc(){}
double TblParamFunc::operator()(ParamSet& X, ParamSet& P){
	return P[0]*m_function(P[1]*X[0]+P[2]);
}