// this file is distributed under 
// GPL v 3.0 license
#ifndef DAIHADDG
#define DAIHADDG
#include <string>
#include <vector>
#include <fit.h>
#define static_right(A) (static_cast<decltype(A)&&>(A))
class hist{
public:
	struct point{
		double x,y,dx,dy;
	};
	hist();
	hist(std::string filename,std::vector<std::string>&&path,std::string histname);
	hist(bool data,std::string reaction,std::vector<std::string>&&path,std::string histname);
	hist(const hist& source);
	hist &operator=(const hist& source);
	virtual ~hist();
	double Entries();
	hist &imbibe(hist& second);
	hist &operator+=(hist& second);
	hist &operator*=(double c);
	hist &operator<<(size_t c);
	hist &operator>>(size_t c);
	point &operator[](int i);
	int count();
	typedef std::vector<point>::iterator iterator;
	typedef std::vector<point>::const_iterator const_iterator;
	iterator begin();
	const_iterator cbegin()const;
	iterator end();
	const_iterator cend() const;
private:
	std::vector<point> data;
	double norm;
};
class PlotHist:public Plot<double>{
public:
    PlotHist();
	PlotHist&Hist(std::string name,hist&&data);
	PlotHist&Hist(std::string name,std::shared_ptr<hist>data);
};
inline std::shared_ptr<Genetic::FitPoints> operator<<(std::shared_ptr<Genetic::FitPoints> dest,hist::point&&source){
	Genetic::FitPoints::Point p;
	p.X<<source.x;p.WX<<source.dx;
	p.y=source.y;p.wy=source.dy;
	return dest<<p;
}
#endif