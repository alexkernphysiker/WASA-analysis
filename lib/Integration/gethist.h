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
		class point{
		public:
			point(const value<double>&pos);
			point(value<double>&&pos);
			point(const value<double>&pos,const value<double>&val);
			point(value<double>&&pos,const value<double>&val);
			point(const value<double>&pos,value<double>&&val);
			point(value<double>&&pos,value<double>&&val);
			point(const point&source);
			value<double>&X()const;
			value<double>&Y()const;
			value<double>&varY();
		private:
			value<double> x;
			value<double> y;
		};
		hist();
		hist(string&&filename,const vector<string>&path,string&&histname);
		hist(string&&filename,vector<string>&&path,string&&histname);
		hist(histsource src,string&&reaction,const vector<string>&path,string&&histname);
		hist(histsource src,string&&reaction,vector<string>&&path,string&&histname);
		hist(const hist& source);
		virtual ~hist();
		hist&operator=(const hist& source);
		
		hist CloneEmptyBins()const;
		double Total()const;
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
		hist &operator+=(const value<double>&c);
		hist &operator+=(value<double>&&c);
		
		hist &operator-=(const hist& second);
		hist &operator-=(function<double(double)>);
		hist &operator-=(const value<double>&c);
		hist &operator-=(value<double>&&c);
		
		hist &operator*=(const hist& second);
		hist &operator*=(function<double(double)>);
		hist &operator*=(const value<double>&c);
		hist &operator*=(value<double>&&c);
		
		hist &operator/=(const hist& second);
		hist &operator/=(function<double(double)>);
		hist &operator/=(const value<double>&c);
		hist &operator/=(value<double>&&c);
		
		hist &operator<<(size_t c);
		hist &operator>>(size_t c);
	private:
		void imbibe(const hist& second);
		vector<point> data;
	};
	inline hist operator+(const hist&a,const hist&b){hist res=a;res+=b;return res;}
	inline hist operator+(hist&&a,const hist&b){hist res=a;res+=b;return res;}
	inline hist operator+(const hist&a,hist&&b){hist res=a;res+=b;return res;}
	inline hist operator+(hist&&a,hist&&b){hist res=a;res+=b;return res;}

	inline hist operator+(const hist&a,function<double(double)>b){hist res=a;res+=b;return res;}
	inline hist operator+(hist&&a,function<double(double)>b){hist res=a;res+=b;return res;}

	inline hist operator+(const hist&a,const value<double>&b){hist res=a;res+=b;return res;}
	inline hist operator+(hist&&a,const value<double>&b){return a+b;}
	inline hist operator+(const hist&a,value<double>&&b){return a+b;}
	inline hist operator+(hist&&a,value<double>&&b){return a+b;}
	
	inline hist operator-(const hist&a,const hist&b){hist res=a;res-=b;return res;}
	inline hist operator-(hist&&a,const hist&b){hist res=a;res-=b;return res;}
	inline hist operator-(const hist&a,hist&&b){hist res=a;res-=b;return res;}
	inline hist operator-(hist&&a,hist&&b){hist res=a;res-=b;return res;}
	
	inline hist operator-(const hist&a,function<double(double)>b){hist res=a;res-=b;return res;}
	inline hist operator-(hist&&a,function<double(double)>b){hist res=a;res-=b;return res;}
	
	inline hist operator-(const hist&a,const value<double>&b){hist res=a;res-=b;return res;}
	inline hist operator-(hist&&a,const value<double>&b){return a-b;}
	inline hist operator-(const hist&a,value<double>&&b){return a-b;}
	inline hist operator-(hist&&a,value<double>&&b){return a-b;}
	
	inline hist operator*(const hist&a,const hist&b){hist res=a;res*=b;return res;}
	inline hist operator*(hist&&a,const hist&b){hist res=a;res*=b;return res;}
	inline hist operator*(const hist&a,hist&&b){hist res=a;res*=b;return res;}
	inline hist operator*(hist&&a,hist&&b){hist res=a;res*=b;return res;}
	
	inline hist operator*(const hist&a,function<double(double)>b){hist res=a;res*=b;return res;}
	inline hist operator*(hist&&a,function<double(double)>b){hist res=a;res*=b;return res;}

	inline hist operator*(const hist&a,const value<double>&b){hist res=a;res*=b;return res;}
	inline hist operator*(hist&&a,const value<double>&b){return a*b;}
	inline hist operator*(const hist&a,value<double>&&b){return a*b;}
	inline hist operator*(hist&&a,value<double>&&b){return a*b;}
	
	inline hist operator/(const hist&a,const hist&b){hist res=a;res/=b;return res;}
	inline hist operator/(hist&&a,const hist&b){hist res=a;res/=b;return res;}
	inline hist operator/(const hist&a,hist&&b){hist res=a;res/=b;return res;}
	inline hist operator/(hist&&a,hist&&b){hist res=a;res/=b;return res;}

	inline hist operator/(const hist&a,function<double(double)>b){hist res=a;res/=b;return res;}
	inline hist operator/(hist&&a,function<double(double)>b){hist res=a;res/=b;return res;}
	
	inline hist operator/(const hist&a,const value<double>&b){hist res=a;res/=b;return res;}
	inline hist operator/(hist&&a,const value<double>&b){return a/b;}
	inline hist operator/(const hist&a,value<double>&&b){return a/b;}
	inline hist operator/(hist&&a,value<double>&&b){return a/b;}
	
	
	double ChiSq(const hist&a,function<double(double)>b,size_t paramcount);
	inline double ChiSq(hist&&a,function<double(double)>b,size_t paramcount){return ChiSq(a,b,paramcount);}
	double ChiSq(const hist&a,const hist&b,size_t paramcount);
	inline double ChiSq(const hist&a,hist&&b,size_t paramcount){return ChiSq(a,b,paramcount);}
	inline double ChiSq(hist&&a,const hist&b,size_t paramcount){return ChiSq(a,b,paramcount);}
	inline double ChiSq(hist&&a,hist&&b,size_t paramcount){return ChiSq(a,b,paramcount);}
	
	class PlotHist:public Plot<double>{
	public:
		PlotHist();
		PlotHist&Hist(const hist&data,string&&title="");
		PlotHist&Hist(hist&&data,string&&title="");
	};
	inline shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints>dest,const hist::point&source){
		return dest<<Point({source.X().val()},{source.X().delta()},source.Y().val(),source.Y().delta());
	}
	inline shared_ptr<FitPoints> operator<<(shared_ptr<FitPoints>dest,const hist&source){
		for(auto&P:source)dest<<P;
		return dest;
	}
};
#endif