// this file is distributed under 
// GPL v 3.0 license
#ifndef FOSfxdVj
#define FOSfxdVj
#include <functional>
#include <abstract.h>
#include <genetic.h>
#include <fit.h>
#include "gethist.h"
template<class GENETIC,class OPTIMALITY,class... ARGS>
class HistFitWithShift:public virtual GENETIC{
public:
	HistFitWithShift(ARGS... args):Genetic::AbstractGenetic(std::make_shared<OPTIMALITY>(args...)),GENETIC(){}
	virtual ~HistFitWithShift(){}
	OPTIMALITY&&GetOpt(){
		return static_cast<OPTIMALITY&&>(*Genetic::AbstractGenetic::OptimalityCalculator());
	}
protected:
	virtual void mutations(Genetic::ParamSet&C,Genetic::RANDOM&R)override{
		double offs=C[0];
		GENETIC::mutations(C,R);
		std::uniform_real_distribution<double> Prob(0,1);
		if(Prob(R)<0.2)
			offs+=1;
		if(Prob(R)<0.2)
			offs-=1;
		C.Set(0,offs);
	}
};
class HistToHists:public Genetic::IOptimalityFunction{
public:
	HistToHists(std::shared_ptr<hist>data,std::vector<std::shared_ptr<hist>>simulation);
    virtual ~HistToHists();
	virtual double operator()(const Genetic::ParamSet&P)const override;
	std::shared_ptr<hist> GetSubHist(size_t i,Genetic::ParamSet&&P);
private:
	std::shared_ptr<hist>m_data;
	std::vector<std::shared_ptr<hist>>m_simulation;
};
#endif 