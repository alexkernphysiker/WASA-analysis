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
	using namespace MathTemplates;
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
	hist<double> ReadHist(const string&filename,const vector<string>&path,const string&histname){
		vector<point<double>> points;
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
					points.push_back(point<double>(value<double>(x,dx),value<double>(y,dy)));
				}
			}
			file->Close();
			delete file;
		}
		return hist<double>(points);
	}
	hist<double> Hist(histsource src, const string&reaction, const vector<string>&path,const string&histname){
		hist<double> res;
		switch(src){
			case MC:{
				res=ReadHist(inputpath+"/MC"+reaction+".root",path,histname);
			}break;
			case DATA:{
				for(ALLRUNS){
					hist<double> tmp=ReadHist(inputpath+"/Data"+reaction+"_run_"+to_string(runindex)+".root",path,histname);
					if(tmp.size()>0)
						if(res.size()==0)
							res=tmp;
						else
							res.imbibe(tmp);
				}
			}break;
		};
		return res;
	}
	hist<double> Hist(histsource src, const string&reaction, const vector<string>&path,string&&histname){return Hist(src,reaction,path,histname);}
	hist<double> Hist(histsource src, string&& reaction, const vector<string>& path, string&& histname){return Hist(src,reaction,path,histname);}
	hist<double> Hist(histsource src, string&& reaction, vector< string >&& path, string&& histname){return Hist(src,reaction,path,histname);}
};