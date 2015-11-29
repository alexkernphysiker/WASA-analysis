// this file is distributed under 
// GPL v 3.0 license
#include <math.h>
#include <TObject.h>
#include <TH1F.h>
#include <TH1.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include "../../phys_constants.h"
#include "gethist.h"
#include "str_get.h"
using namespace std;
using namespace Genetic;
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
hist::hist(string&&filename,vector<string>&&path,string&&histname){
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
hist::hist(histsource src,string&&reaction,vector<string>&&path,string&&histname){
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
			throw math_h_error<hist>("Cannot imbibe histogram. bins differ");
	}
}
hist& hist::operator+=(const hist& second){
	for(size_t i=0,n=count();i<n;i++){
		if(data[i].x==second[i].x){
			data[i].y+=second[i].y;
			data[i].dy+=second[i].dy;
		}else
			throw math_h_error<hist>("Cannot add histogram. bins differ");
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
			data[i].dy+=second[i].dy;
		}else
			throw math_h_error<hist>("Cannot substract histogram. bins differ");
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
			data[i].dy=Y1.second*Y2.first+Y2.second*Y1.first;
		}else
			throw math_h_error<hist>("Cannot multiply by a histogram. bins differ");
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
			data[i].dy=Y1.second/Y2.first+Y2.second*Y1.first/pow(Y2.first,2);
		}else
			throw math_h_error<hist>("Cannot multiply by a histogram. bins differ");
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
double hist::HowClose(const hist& second_hist) const{
	double res=0;
	for(size_t i=1,n=count()-1;i<n;i++){
		if(data[i].x==second_hist[i].x){
			auto pnt=[i](const hist&H){
				return make_pair(
					H.data[i].y,
					sqrt(pow(H.data[i].dy,2)+pow(H.data[i+1].y-H.data[i-1].y,2)/4.0)
				);
			};
			auto Y1=pnt(*this),Y2=pnt(second_hist);
			res+=pow(Y1.first-Y2.first,2)/(pow(Y1.second,2)+pow(Y2.second,2));
		}else
			throw math_h_error<hist>("Cannot calculate optimality of fit by another histogram. bins differ");
	}
	return res;
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

