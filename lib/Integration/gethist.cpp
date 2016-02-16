// this file is distributed under 
// MIT license
#include <math.h>
#include <TObject.h>
#include <TH1F.h>
#include <TH1.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <math_h/error.h>
#include "experiment_conv.h"
#include "gethist.h"
#include "str_get.h"
namespace ROOT_data{
	using namespace std;
	using namespace Genetic;
	using namespace MathTemplates;
	using namespace GnuplotWrap;
	string inputpath=ENV(PRESEL_DATA);
	string outpath=ENV(OUTPUT_PLOTS);
	double PresentRunsAmountRatio(string&&reaction){
		size_t allruns=0,present_runs=0;
		for(ALLRUNS){
			allruns++;
			TFile *file=TFile::Open((inputpath+"/Data"+reaction+"_run_"+to_string(runindex)+".root").c_str());
			if(file){
				present_runs++;
				file->Close();
				delete file;
			}
		}
		return double(present_runs)/double(allruns);
	}
	hist::point::point(const value<double>& pos):x(pos){}
	hist::point::point(value<double>&& pos):x(pos){}
	hist::point::point(const value<double>& pos, const value<double>& val):x(pos),y(val){}
	hist::point::point(value<double>&& pos, const value<double>& val):x(pos),y(val){}
	hist::point::point(const value<double>& pos, value<double>&& val):x(pos),y(val){}
	hist::point::point(value<double>&& pos, value<double>&& val):x(pos),y(val){}
	hist::point::point(const hist::point& source):x(source.x),y(source.y){}
	value<double>& hist::point::X() const{return const_cast<value<double>&>(x);}
	value<double>& hist::point::Y() const{return const_cast<value<double>&>(y);}
	value<double>& hist::point::varY(){return y;}
	hist::hist(){}
	hist::hist(string&&filename,const vector<string>&path,string&&histname){
		TFile* file=TFile::Open(filename.c_str());
		if(file){
			TDirectoryFile* dir1=file;
			for(string name:path){
				TDirectoryFile* dir2=dynamic_cast<TDirectoryFile*>(dir1->Get(name.c_str()));
				if(dir2)
					dir1=dir2;
			}
			TH1F* histogram=dynamic_cast<TH1F*>(dir1->Get(histname.c_str()));
			if(histogram){
				for(int i=1,N=histogram->GetNbinsX();i<=N;i++){
					double y=histogram->GetBinContent(i);
					double dy=sqrt(y);
					if(dy<1.0)
						dy=1.0;
					double x=histogram->GetBinCenter(i);
					double dx=histogram->GetBinWidth(i)/2.0;
					data.push_back(point(value<double>(x,dx),value<double>(y,dy)));
				}
			}
			file->Close();
			delete file;
		}
	}
	hist::hist(string&& filename, vector< string >&& path, string&& histname):hist(static_cast<string&&>(filename),path,static_cast<string&&>(histname)){}
	hist::hist(histsource src, const string& reaction, const vector< string >& path, string&& histname){
		switch(src){
			case MC:{
				hist tmp(inputpath+"/MC"+reaction+".root",static_cast<decltype(path)&&>(path),static_cast<string&&>(histname));
				operator=(tmp);
			}break;
			case DATA:{
				data.clear();
				for(ALLRUNS){
					hist tmp(inputpath+"/Data"+reaction+"_run_"+to_string(runindex)+".root",static_cast<decltype(path)&&>(path),static_cast<string&&>(histname));
					if(tmp.data.size()>0){
						if(data.size()==0)
							operator=(tmp);
						else
							imbibe(tmp);
					}
				}
			}break;
		};
	}
	hist::hist(histsource src, string&& reaction, const vector< string >& path, string&& histname):hist(src,reaction,path,static_cast<string&&>(histname)){}
	hist::hist(histsource src,string&&reaction,vector<string>&&path,string&&histname):hist(src,reaction,path,static_cast<string&&>(histname)){}
	hist::hist(const hist& source){
		for(point p: source.data)
			data.push_back(p);
	}
	hist& hist::Cut(double a, double b){
		for(int i=0;i<data.size();i++)if((data[i].X().val()<a)||(data[i].X().val()>b)){
			data.erase(data.begin()+i);
			i--;
		}
		return *this;
	}
	
	
	hist hist::CloneEmptyBins()const{
		hist res(*this);
		for(int i=0,n=count();i<n;i++){
			res.data[i].varY()=value<double>(0,0);
		}
		return res;
	}
	
	hist& hist::operator=(const hist& source){
		data.clear();
		for(point p: source)
			data.push_back(p);
		return *this;
	}
	hist::~hist(){}
	size_t hist::count()const{
		return data.size();
	}
	double hist::Total()const{
		double res=0;
		for(const point& P:data)
			res+=P.Y().val();
		return res;
	}
	hist::point& hist::operator[](int i)const{return const_cast<point&>(data[i]);}
	hist::point& hist::operator()(double x)const{
		size_t res_ind=0;
		double difference=+INFINITY;
		for(size_t i=0,n=count();i<n;i++){
			double d=pow(data[i].X().val()-x,2);
			if(d<difference){
				difference=d;
				res_ind=i;
			}
		}
		return operator[](res_ind);
	}
	
	hist::iterator hist::begin(){return data.begin();}
	hist::const_iterator hist::begin() const{return data.begin();}
	hist::const_iterator hist::cbegin() const{return data.cbegin();}
	hist::iterator hist::end(){return data.end();}
	hist::const_iterator hist::end() const{return data.end();}
	hist::const_iterator hist::cend() const{return data.cend();}
	
	void hist::imbibe(const hist& second){
		for(int i=0,n=count();i<n;i++){
			if(data[i].X().val()==second[i].X().val()){
				data[i].varY()=(data[i].Y().val()+second[i].Y().val());
			}else
				throw Exception<hist>("Cannot imbibe histogram. bins differ");
		}
	}
	hist& hist::operator+=(const hist& second){
		for(size_t i=0,n=count();i<n;i++){
			if(data[i].X().val()==second[i].X().val()){
				data[i].varY()+=second[i].Y();
			}else
				throw Exception<hist>("Cannot add histogram. bins differ");
		}
		return *this;
	}
	hist& hist::operator+=(function<double(double)>f){
		for(size_t i=0,n=count();i<n;i++)
			data[i].varY()+=value<double>(f(data[i].X().val()),0);
		return *this;
	}
	hist& hist::operator+=(const value<double>&v){
		for(size_t i=0,n=count();i<n;i++)data[i].varY()+=v;
		return *this;
	}
	hist& hist::operator+=(value<double>&&v){return operator+=(v);}
	
	hist& hist::operator-=(const hist& second){
		for(size_t i=0,n=count();i<n;i++){
			if(data[i].X().val()==second[i].X().val()){
				data[i].varY()-=second[i].Y();
			}else
				throw Exception<hist>("Cannot substract histogram. bins differ");
		}
		return *this;
	}
	hist& hist::operator-=(function<double(double)>f){
		for(size_t i=0,n=count();i<n;i++)
			data[i].varY()-=value<double>(f(data[i].X().val()),0);
		return *this;
	}
	hist& hist::operator-=(const value<double>&v){
		for(size_t i=0,n=count();i<n;i++)data[i].varY()-=v;
		return *this;
	}
	hist& hist::operator-=(value<double>&&v){return operator-=(v);}
	
	
	hist& hist::operator*=(const hist& second){
		for(size_t i=0,n=count();i<n;i++){
			if(data[i].X().val()==second[i].X().val()){
				data[i].varY()*=second[i].Y();
			}else
				throw Exception<hist>("Cannot substract histogram. bins differ");
		}
		return *this;
	}
	hist& hist::operator*=(function<double(double)>f){
		for(size_t i=0,n=count();i<n;i++)
			data[i].varY()*=value<double>(f(data[i].X().val()),0);
		return *this;
	}
	hist& hist::operator*=(const value<double>&v){
		for(size_t i=0,n=count();i<n;i++)data[i].varY()*=v;
		return *this;
	}
	hist& hist::operator*=(value<double>&&v){return operator*=(v);}
	
	
	hist& hist::operator/=(const hist& second){
		for(size_t i=0,n=count();i<n;i++){
			if(data[i].X().val()==second[i].X().val()){
				data[i].varY()/=second[i].Y();
			}else
				throw Exception<hist>("Cannot substract histogram. bins differ");
		}
		return *this;
	}
	hist& hist::operator/=(function<double(double)>f){
		for(size_t i=0,n=count();i<n;i++)
			data[i].varY()/=value<double>(f(data[i].X().val()),0);
		return *this;
	}
	hist& hist::operator/=(const value<double>&v){
		for(size_t i=0,n=count();i<n;i++)data[i].varY()/=v;
		return *this;
	}
	hist& hist::operator/=(value<double>&&v){return operator/=(v);}
	
	
	
	hist& hist::operator<<(size_t c){
		for(size_t i=0,n=count()-c;i<n;i++)
			data[i].varY()=data[i+c].Y();
		for(size_t n=count(),i=n-c;i<n;i++)
			data[i].varY()=value<double>(0);
		return *this;
	}
	hist& hist::operator>>(size_t c){
		for(size_t i=count()-1;i>=c;i--)
			data[i].varY()=data[i+c].Y();
		for(size_t i=0;i<c;i++)
			data[i].varY()=value<double>(0);
		return *this;
	}
	
	double ChiSq(const hist&a,function<double(double)>b,size_t paramcount){
		double res=0,k=a.count()-paramcount;
		if(k<=0)throw Exception<hist>("ChiSq error: too few points or too many parameters");
		for(const hist::point&p:a)
			res+=pow((p.Y().val()-b(p.X().val()))/p.Y().val(),2);
		return res/k;
	}
	double ChiSq(const hist&a,const hist&b,size_t paramcount){
		if(a.count()!=b.count())throw Exception<hist>("ChiSq error: hists size mismatch");
		double res=0,k=a.count()-paramcount;
		if(k<=0)throw Exception<hist>("ChiSq error: too few points or too many parameters");
		for(size_t i=0,n=a.count();i<n;i++){
			if(a[i].X().val()!=b[i].X().val())throw Exception<hist>("ChiSq error: hists bins mismatch");
			if(a[i].X().delta()!=b[i].X().delta())throw Exception<hist>("ChiSq error: hists bins widths mismatch");
			auto diff=a[i].Y()-b[i].Y();
			res+=pow(diff.val()/diff.delta(),2);
		}
		return res/k;
	}
	PlotHist::PlotHist():Plot<double>(){}
	PlotHist& PlotHist::Hist(const hist&data,const string&title){
		Plot<double>::OutputPlot([&data](std::ofstream&str){
			for(hist::point p:data)
				str<<p.X().val()<<" "<<p.Y().val()<<" "<<p.X().delta()<<" "<<p.Y().delta()<<endl;
		},"using 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars",title);
		return *this;
	}
	PlotHist& PlotHist::Hist(const hist& data, string&& title){return Hist(data,title);}
	PlotHist& PlotHist::Hist(hist&& data,string&&title){return Hist(data,title);}
};