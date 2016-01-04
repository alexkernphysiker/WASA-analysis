// this file is distributed under 
// MIT license
#ifndef DAIHADDG
#define DAIHADDG
#include <string>
#include <functional>
#include <vector>
#include <Genetic/fit.h>
namespace ROOT_data{
	using namespace std;
	using namespace Genetic;
	double PresentRunsAmountRatio(string&&reaction);
	enum histsource{MC,DATA};
	class hist{
	public:
		struct point{
			double x,y,dx,dy;
		};
		hist();
		hist(string&&filename,vector<string>&&path,string&&histname);
		hist(histsource src,string&&reaction,vector<string>&&path,string&&histname);
		hist(const hist& source);
		virtual ~hist();
		hist&operator=(const hist& source);
		
		hist CloneEmptyBins()const;
		double Entries()const;
		point&operator[](int i)const;
		point&operator()(double x)const;
		size_t count()const;
		
		typedef vector<point>::iterator iterator;
		typedef vector<point>::const_iterator const_iterator;
		iterator begin();
		const_iterator begin()const;
		const_iterator cbegin()const;
		iterator end();
		const_iterator end() const;
		const_iterator cend() const;
		
		hist& Cut(double a,double b);
		hist &operator+=(const hist& second);
		hist &operator+=(function<double(double)>);
		hist &operator-=(const hist& second);
		hist &operator-=(function<double(double)>);
		hist &operator*=(const hist& second);
		hist &operator*=(double c);
		hist &operator*=(function<double(double)>);
		hist &operator/=(const hist& second);
		hist &operator/=(double c);
		hist &operator/=(function<double(double)>);
		hist &operator<<(size_t c);
		hist &operator>>(size_t c);
	private:
		void imbibe(const hist& second);
		vector<point> data;
	};
	class PlotHist:public Plot<double>{
	public:
		PlotHist();
		PlotHist&Hist(string&&name,const hist&data);
		PlotHist&HistWLine(string&&name,const hist&data);
	};
	inline shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints>dest,const hist::point&source){
		return dest<<Point({source.x},{source.dx},source.y,source.dy);
	}
	inline shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints>dest,const hist&source){
		for(auto&P:source)dest<<P;
		return dest;
	}
};
#endif