#ifndef WJUPOASA
#define WJUPOASA
#include <paramfunc.h>
#include <math_h/interpolate.h>
typedef LinearInterpolation<double> FuncTbl;
class TblParamFunc:public virtual Genetic::IParamFunc{
public:
	TblParamFunc(FuncTbl func);
	virtual ~TblParamFunc();
	virtual double operator()(Genetic::ParamSet& X,Genetic::ParamSet& P) override;
	enum{ParamCount=3,ArgCount=1};
private:
	FuncTbl m_function;
};

#endif