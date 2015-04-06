#include <math.h>
#include <TObject.h>
#include <TH1F.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include "gethist.h"
using namespace std;
hist::hist(std::string filename,std::vector<std::string> &path,std::string histname){
	TFile* file=TFile::Open(filename.c_str());
	if(file){
		TDirectoryFile* dir1=file;
		for(string name:path){
			TDirectoryFile* dir2=dynamic_cast<TDirectoryFile*>(dir1->Get(name.c_str()));
			if(dir2)
				dir1=dir2;
			else
				throw;
		}
		TH1F* histogram=dynamic_cast<TH1F*>(dir1->Get(histname.c_str()));
		if(histogram){
			norm=histogram->GetEntries();
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
		}else
			throw;
		
	}else
		throw;
}
hist::~hist(){}
int hist::count(){
	return data.size();
}
double hist::Entries(){
	return norm;
}
point& hist::operator[](int i){
	return data[i];
}
hist::iterator hist::begin(){
	return data.begin();
}
hist::const_iterator hist::cbegin() const{
	return data.cbegin();
}
hist::iterator hist::end(){
	return data.end();
}
hist::const_iterator hist::cend() const{
	return data.cend();
}


