// this file is distributed under 
// GPL v 3.0 license
#ifndef DAIHADDG
#define DAIHADDG
#include <string>
#include <functional>
#include <vector>
#include <fit.h>
double PresentRunsAmountRatio(std::string&&reaction);
class hist{
public:
	struct point{
		double x,y,dx,dy;
	};
	hist();
	hist(std::string&&filename,std::vector<std::string>&&path,std::string&&histname);
	hist(bool data,std::string&&reaction,std::vector<std::string>&&path,std::string&&histname);
	hist(const hist& source);
	virtual ~hist();
	hist& Cut(double a,double b);

	hist&operator=(const hist& source);
	double HowClose(const hist&B)const;
	hist CloneEmptyBins()const;
	double Entries()const;
	point&operator[](int i)const;
	int count()const;
	typedef std::vector<point>::iterator iterator;
	typedef std::vector<point>::const_iterator const_iterator;
	iterator begin();
	const_iterator begin()const;
	const_iterator cbegin()const;
	iterator end();
	const_iterator end() const;
	const_iterator cend() const;
	
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
	double norm;
};
class PlotHist:public Plot<double>{
public:
    PlotHist();
	PlotHist& Hist(std::string&&name,const hist&data);
	PlotHist& HistWLine(std::string&&name,const hist&data);
};
inline std::shared_ptr<Genetic::FitPoints> operator<<(std::shared_ptr<Genetic::FitPoints> dest,const hist::point&source){
	Genetic::FitPoints::Point p;
	p.X<<source.x;p.WX<<source.dx;
	p.y=source.y;p.wy=source.dy;
	return dest<<p;
}
#endif