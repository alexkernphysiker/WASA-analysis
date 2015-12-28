// this file is distributed under 
// MIT license
#ifndef DAIHADDG
#define DAIHADDG
#include <string>
#include <functional>
#include <vector>
#include <Genetic/fit.h>
double PresentRunsAmountRatio(std::string&&reaction);
enum histsource{MC,DATA};
class hist{
public:
	struct point{
		double x,y,dx,dy;
		
	};
	hist();
	hist(std::string&&filename,std::vector<std::string>&&path,std::string&&histname);
	hist(histsource src,std::string&&reaction,std::vector<std::string>&&path,std::string&&histname);
	hist(const hist& source);
	virtual ~hist();
	hist&operator=(const hist& source);

	hist CloneEmptyBins()const;
	double Entries()const;
	point&operator[](int i)const;
	point&operator()(double x)const;
	size_t count()const;

	typedef std::vector<point>::iterator iterator;
	typedef std::vector<point>::const_iterator const_iterator;
	iterator begin();
	const_iterator begin()const;
	const_iterator cbegin()const;
	iterator end();
	const_iterator end() const;
	const_iterator cend() const;
	
	hist& Cut(double a,double b);
	hist &operator+=(const hist& second);
	hist &operator+=(std::function<double(double)>);
	hist &operator-=(const hist& second);
	hist &operator-=(std::function<double(double)>);
	hist &operator*=(const hist& second);
	hist &operator*=(double c);
	hist &operator*=(std::function<double(double)>);
	hist &operator/=(const hist& second);
	hist &operator/=(double c);
	hist &operator/=(std::function<double(double)>);
	hist &operator<<(size_t c);
	hist &operator>>(size_t c);
private:
	void imbibe(const hist& second);
	std::vector<point> data;
};
class PlotHist:public Plot<double>{
public:
    PlotHist();
	PlotHist& Hist(std::string&&name,const hist&data);
	PlotHist& HistWLine(std::string&&name,const hist&data);
};
inline std::shared_ptr<Genetic::FitPoints> operator<<(std::shared_ptr<Genetic::FitPoints>dest,const hist::point&source){
	return dest<<Genetic::Point({source.x},{source.dx},source.y,source.dy);
}
inline std::shared_ptr<Genetic::FitPoints> operator<<(std::shared_ptr<Genetic::FitPoints>dest,const hist&source){
	for(auto&P:source)
		dest<<P;
	return dest;
}
#endif