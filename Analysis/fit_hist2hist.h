// this file is distributed under 
// GPL v 3.0 license
#ifndef FOSfxdVj
#define FOSfxdVj
#include <functional>
#include <abstract.h>
#include "gethist.h"
class Hist2Hist:public Genetic::IOptimalityFunction{
public:
	typedef std::function<double(double,Genetic::ParamSet&&)> bg_func;
	Hist2Hist(std::shared_ptr<hist>data,std::shared_ptr<hist>mc,bg_func background);
    virtual ~Hist2Hist();
    virtual double operator()(Genetic::ParamSet&& P)final;
private:
	std::shared_ptr<hist> DATA,MC;
	bg_func BG;
};
template<class GENETIC>
class FitHist:public virtual GENETIC{
public:
	FitHist(std::shared_ptr<hist>data,std::shared_ptr<hist>mc,Hist2Hist::bg_func background)
		:Genetic::AbstractGenetic(std::make_shared<Hist2Hist>(data,mc,background)),GENETIC(){}
};
#endif 