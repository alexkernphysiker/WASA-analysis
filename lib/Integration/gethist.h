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
	class value{
	public:
		value();
		value(double v);
		value(double v,double err);
		value(const value&source);
		double val()const;
		double delta()const;
		double epsilon()const;
		double min()const;
		double max()const;
		value&operator+=(const value&other);
		value&operator-=(const value&other);
		value&operator*=(const value&other);
		value&operator/=(const value&other);
	private:
		double Value,Error;
	};
	inline value operator+(const value&a,const value&b){value res=a;res+=b;return res;}
	inline value operator+(value&&a,const value&b){return a+b;}
	inline value operator+(const value&a,value&&b){return a+b;}
	inline value operator+(value&&a,value&&b){return a+b;}
	inline value operator-(const value&a,const value&b){value res=a;res-=b;return res;}
	inline value operator-(value&&a,const value&b){return a-b;}
	inline value operator-(const value&a,value&&b){return a-b;}
	inline value operator-(value&&a,value&&b){return a-b;}
	inline value operator*(const value&a,const value&b){value res=a;res*=b;return res;}
	inline value operator*(value&&a,const value&b){return a*b;}
	inline value operator*(const value&a,value&&b){return a*b;}
	inline value operator*(value&&a,value&&b){return a*b;}
	inline value operator/(const value&a,const value&b){value res=a;res/=b;return res;}
	inline value operator/(value&&a,const value&b){return a/b;}
	inline value operator/(const value&a,value&&b){return a/b;}
	inline value operator/(value&&a,value&&b){return a/b;}
	class hist{
	public:
		class point{
		public:
			point(const value&pos);
			point(value&&pos);
			point(const value&pos,const value&val);
			point(value&&pos,const value&val);
			point(const value&pos,value&&val);
			point(value&&pos,value&&val);
			point(const point&source);
			value&X()const;
			value&Y()const;
			value&varY();
		private:
			value x;
			value y;
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
		hist &operator+=(const value&c);
		hist &operator+=(value&&c);
		
		hist &operator-=(const hist& second);
		hist &operator-=(function<double(double)>);
		hist &operator-=(const value&c);
		hist &operator-=(value&&c);
		
		hist &operator*=(const hist& second);
		hist &operator*=(function<double(double)>);
		hist &operator*=(const value&c);
		hist &operator*=(value&&c);
		
		hist &operator/=(const hist& second);
		hist &operator/=(function<double(double)>);
		hist &operator/=(const value&c);
		hist &operator/=(value&&c);
		
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

	inline hist operator+(const hist&a,const value&b){hist res=a;res+=b;return res;}
	inline hist operator+(hist&&a,const value&b){return a+b;}
	inline hist operator+(const hist&a,value&&b){return a+b;}
	inline hist operator+(hist&&a,value&&b){return a+b;}
	
	inline hist operator-(const hist&a,const hist&b){hist res=a;res-=b;return res;}
	inline hist operator-(hist&&a,const hist&b){hist res=a;res-=b;return res;}
	inline hist operator-(const hist&a,hist&&b){hist res=a;res-=b;return res;}
	inline hist operator-(hist&&a,hist&&b){hist res=a;res-=b;return res;}
	
	inline hist operator-(const hist&a,function<double(double)>b){hist res=a;res-=b;return res;}
	inline hist operator-(hist&&a,function<double(double)>b){hist res=a;res-=b;return res;}
	
	inline hist operator-(const hist&a,const value&b){hist res=a;res-=b;return res;}
	inline hist operator-(hist&&a,const value&b){return a-b;}
	inline hist operator-(const hist&a,value&&b){return a-b;}
	inline hist operator-(hist&&a,value&&b){return a-b;}
	
	inline hist operator*(const hist&a,const hist&b){hist res=a;res*=b;return res;}
	inline hist operator*(hist&&a,const hist&b){hist res=a;res*=b;return res;}
	inline hist operator*(const hist&a,hist&&b){hist res=a;res*=b;return res;}
	inline hist operator*(hist&&a,hist&&b){hist res=a;res*=b;return res;}
	
	inline hist operator*(const hist&a,function<double(double)>b){hist res=a;res*=b;return res;}
	inline hist operator*(hist&&a,function<double(double)>b){hist res=a;res*=b;return res;}

	inline hist operator*(const hist&a,const value&b){hist res=a;res*=b;return res;}
	inline hist operator*(hist&&a,const value&b){return a*b;}
	inline hist operator*(const hist&a,value&&b){return a*b;}
	inline hist operator*(hist&&a,value&&b){return a*b;}
	
	inline hist operator/(const hist&a,const hist&b){hist res=a;res/=b;return res;}
	inline hist operator/(hist&&a,const hist&b){hist res=a;res/=b;return res;}
	inline hist operator/(const hist&a,hist&&b){hist res=a;res/=b;return res;}
	inline hist operator/(hist&&a,hist&&b){hist res=a;res/=b;return res;}

	inline hist operator/(const hist&a,function<double(double)>b){hist res=a;res/=b;return res;}
	inline hist operator/(hist&&a,function<double(double)>b){hist res=a;res/=b;return res;}
	
	inline hist operator/(const hist&a,const value&b){hist res=a;res/=b;return res;}
	inline hist operator/(hist&&a,const value&b){return a/b;}
	inline hist operator/(const hist&a,value&&b){return a/b;}
	inline hist operator/(hist&&a,value&&b){return a/b;}
	
	
	double ChiSq(const hist&a,function<double(double)>b,size_t paramcount);
	inline double ChiSq(hist&&a,function<double(double)>b,size_t paramcount){return ChiSq(a,b,paramcount);}
	double ChiSq(const hist&a,const hist&b,size_t paramcount);
	inline double ChiSq(const hist&a,hist&&b,size_t paramcount){return ChiSq(a,b,paramcount);}
	inline double ChiSq(hist&&a,const hist&b,size_t paramcount){return ChiSq(a,b,paramcount);}
	inline double ChiSq(hist&&a,hist&&b,size_t paramcount){return ChiSq(a,b,paramcount);}
	
	class PlotHist:public Plot<double>{
	public:
		PlotHist();
		PlotHist&Hist(string&&name,const hist&data);
		PlotHist&Hist(string&&name,hist&&data);
		PlotHist&HistWLine(string&&name,const hist&data);
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