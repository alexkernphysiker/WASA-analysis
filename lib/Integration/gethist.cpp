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
#include "../../phys_constants.h"
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
					point p;
					p.x=x;p.y=y;p.dx=dx;p.dy=dy;
					data.push_back(p);
				}
			}
			file->Close();
			delete file;
		}
	}
	hist::hist(string&& filename, vector< string >&& path, string&& histname):hist(static_cast<string&&>(filename),path,static_cast<string&&>(histname)){}
	hist::hist(histsource src, string&& reaction, const vector< string >& path, string&& histname){
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
	hist::hist(histsource src,string&&reaction,vector<string>&&path,string&&histname):hist(src,static_cast<string&&>(reaction),path,static_cast<string&&>(histname)){}
	hist::hist(const hist& source){
		for(point p: source.data)
			data.push_back(p);
	}
	hist& hist::Cut(double a, double b){
		for(int i=0;i<data.size();i++)if((data[i].x<a)||(data[i].x>b)){
			data.erase(data.begin()+i);
			i--;
		}
		return *this;
	}
	
	
	hist hist::CloneEmptyBins()const{
		hist res(*this);
		for(int i=0,n=count();i<n;i++){
			res.data[i].y=0;
			res.data[i].dy=1;
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
	double hist::Entries()const{
		double res=0;
		for(const point& P:data)
			res+=P.y;
		return res;
	}
	hist::point& hist::operator[](int i)const{return const_cast<point&>(data[i]);}
	hist::point& hist::operator()(double x)const{
		size_t res_ind=0;
		double difference=+INFINITY;
		for(size_t i=0,n=count();i<n;i++){
			double d=pow(data[i].x-x,2);
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
			if(data[i].x==second[i].x){
				data[i].y+=second[i].y;
				data[i].dy=sqrt(data[i].y);
				if(data[i].dy<1)data[i].dy=1;
			}else
				throw Exception<hist>("Cannot imbibe histogram. bins differ");
		}
	}
	hist& hist::operator+=(const hist& second){
		for(size_t i=0,n=count();i<n;i++){
			if(data[i].x==second[i].x){
				data[i].y+=second[i].y;
				data[i].dy=sqrt(pow(data[i].dy,2)+pow(second[i].dy,2));
			}else
				throw Exception<hist>("Cannot add histogram. bins differ");
		}
		return *this;
	}
	hist& hist::operator+=(function<double(double)>f){
		for(size_t i=0,n=count();i<n;i++)
			data[i].y+=f(data[i].x);
		return *this;
	}
	hist& hist::operator-=(const hist& second){
		for(size_t i=0,n=count();i<n;i++){
			if(data[i].x==second[i].x){
				data[i].y-=second[i].y;
				data[i].dy=sqrt(pow(data[i].dy,2)+pow(second[i].dy,2));
			}else
				throw Exception<hist>("Cannot substract histogram. bins differ");
		}
		return *this;
	}
	hist& hist::operator-=(function<double(double)>f){
		for(size_t i=0,n=count();i<n;i++)
			data[i].y-=f(data[i].x);
		return *this;
	}
	hist& hist::operator*=(const hist& second_hist){
		for(size_t i=0,n=count();i<n;i++){
			if(data[i].x==second_hist[i].x){
				auto Y1=make_pair(data[i].y,data[i].dy);
				auto Y2=make_pair(second_hist[i].y,second_hist[i].dy);
				data[i].y=Y1.first*Y2.first;
				data[i].dy=sqrt(pow(Y1.second*Y2.first,2)+pow(Y2.second*Y1.first,2));
			}else
				throw Exception<hist>("Cannot multiply by a histogram. bins differ");
		}
		return *this;
	}
	hist& hist::operator*=(double c){
		for(size_t i=0,n=count();i<n;i++){
			data[i].y*=c;
			if(c>=0)data[i].dy*=c;
			else data[i].dy*=-c;
		}
		return *this;
	}
	hist& hist::operator*=(function<double(double)>f){
		for(int i=0,n=count();i<n;i++){
			double c=f(data[i].x);
			data[i].y*=c;
			if(c>=0)data[i].dy*=c;
			else data[i].dy*=-c;
		}
		return *this;
	}
	hist& hist::operator/=(const hist& second_hist){
		for(size_t i=0,n=count();i<n;i++){
			if(data[i].x==second_hist[i].x){
				auto Y1=make_pair(data[i].y,data[i].dy);
				auto Y2=make_pair(second_hist[i].y,second_hist[i].dy);
				data[i].y=Y1.first/Y2.first;
				data[i].dy=sqrt(pow(Y1.second/Y2.first,2)+pow(Y2.second*Y1.first/pow(Y2.first,2),2));
			}else
				throw Exception<hist>("Cannot multiply by a histogram. bins differ");
		}
		return *this;
	}
	hist& hist::operator/=(double c){
		return operator*=(1.0/c);
	}
	hist& hist::operator/=(function<double(double)>f){
		return operator*=([f](double x){return 1.0/f(x);});
	}
	hist& hist::operator<<(size_t c){
		for(size_t i=0,n=count()-c;i<n;i++){
			data[i].y=data[i+c].y;
			data[i].dy=data[i+c].dy;
		}
		for(size_t n=count(),i=n-c;i<n;i++){
			data[i].y=0;
			data[i].dy=1;
		}
		return *this;
	}
	hist& hist::operator>>(size_t c){
		for(size_t i=count()-1;i>=c;i--){
			data[i].y=data[i-c].y;
			data[i].dy=data[i-c].dy;
		}
		for(size_t i=0;i<c;i++){
			data[i].y=0;
			data[i].dy=1;
		}
		return *this;
	}
	
	double ChiSq(const hist&a,function<double(double)>b,size_t paramcount){
		double res=0,k=a.count()-paramcount;
		if(k<=0)throw Exception<hist>("ChiSq error: too few points or too many parameters");
		for(const hist::point&p:a)
			res+=pow((p.y-b(p.x))/p.dy,2);
		return res/k;
	}
	double ChiSq(const hist&a,const hist&b,size_t paramcount){
		if(a.count()!=b.count())throw Exception<hist>("ChiSq error: hists size mismatch");
		double res=0,k=a.count()-paramcount;
		if(k<=0)throw Exception<hist>("ChiSq error: too few points or too many parameters");
		for(size_t i=0,n=a.count();i<n;i++){
			if(a[i].x!=b[i].x)throw Exception<hist>("ChiSq error: hists bins mismatch");
			if(a[i].dx!=b[i].dx)throw Exception<hist>("ChiSq error: hists bins widths mismatch");
			res+=pow(a[i].y-b[i].x,2)/(pow(a[i].dy,2)+pow(b[i].dy,2));
		}
		return res/k;
	}
	
	PlotHist::PlotHist():Plot<double>(){}
	PlotHist& PlotHist::Hist(string&&name,const hist&data){
		Plot<double>::OutputPlot(static_cast<string&&>(name),[&data](std::ofstream&str){
			for(hist::point p:data)str<<p.x<<" "<<p.y<<" "<<p.dx<<" "<<p.dy<<"\n";
		},"using 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars");
		return *this;
	}
	PlotHist& PlotHist::HistWLine(string&& name, const hist& data){
		Plot<double>::OutputPlot(static_cast<string&&>(name),[&data](std::ofstream&str){
			for(hist::point p:data)str<<p.x<<" "<<p.y<<" "<<p.dx<<" "<<p.dy<<"\n";
		},"using 1:2 w l");
		return *this;
	}
};