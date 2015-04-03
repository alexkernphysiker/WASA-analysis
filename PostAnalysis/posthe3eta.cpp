#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <memory>
#include <unistd.h>

#include <TObject.h>
#include <TH1F.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectoryFile.h>

#include <paramfunc.h>
#include <fitpoints.h>
#include <filter.h>
#include <initialconditions.h>
#include <genetic.h>

#include <phys_constants.h>
using namespace std;
using namespace Fit;
int main(int arg_cnt,char**arg) {
	printf("This will perform the post analysis of data from p d -> He3 eta reaction\n");
	string inputpath;{
		stringbuf buffer;
		ostream os (&buffer); 
		os<<getenv("WASA_OUTPUT_DATA");
		inputpath=buffer.str();
		printf("input path: %s\n",inputpath.c_str());
	}
	string outpath;{
		stringbuf buffer;
		ostream os (&buffer); 
		os<<getenv("POST_ANALYSIS_DATA");
		outpath=buffer.str();
		printf("output path: %s\n",outpath.c_str());
	}
	auto points_mc=make_shared<ChiSquare>();{
		string MCFile=inputpath+"/MCHe3Eta.root";
		TFile* file=TFile::Open(MCFile.c_str());
		if(file){
			TDirectoryFile* dir1=dynamic_cast<TDirectoryFile*>(file->Get("Histograms"));
			if(dir1){
				TDirectoryFile* dir2=dynamic_cast<TDirectoryFile*>(dir1->Get("MissingMass"));
				if(dir2){
					TH1F* histogram=dynamic_cast<TH1F*>(dir2->Get("MissingMass"));
					if(histogram){
						double norm=histogram->GetEntries();
						printf("Hist: %f entries\n",norm);
						for(int i=1,N=histogram->GetNbinsX();i<=N;i++){
							double y=histogram->GetBinContent(i);
							double dy=sqrt(y);
							if(dy<1)dy=1;
							double x=histogram->GetBinCenter(i);
							double dx=histogram->GetBinWidth(i)/2.0;
							points_mc->Add(ParamSet(x),ParamSet(dx),y/norm,dy/norm);
						}
					}else
						printf("No hist\n");
				}else
					printf("No tree2\n");
			}else
				printf("No tree1\n");
		}else
			printf("No file\n");
	}
	printf("Fitting MonteCarlo...\n");
	points_mc=SelectFitPoints(points_mc,make_shared<Filter<>>([](ParamSet&X){
		return (X[0]>0.52)&&(X[0]<0.57);
	}));
	DifferentialRandomMutations<> fit_mc(
		make_shared<Mul<Func3<Gaussian,Arg<0>,Par<0>,Par<1>>,Par<2>>>(),
		points_mc,THREADS_COUNT
	);
	fit_mc.SetFilter(make_shared<And>()
		<<(make_shared<Above>()<<0.53<<0<<0)
		<<(make_shared<Below>()<<0.56)
	);
	fit_mc.Init(60,make_shared<GenerateByGauss>()
		<<make_pair(m_eta,0.05)<<make_pair(0.05,0.05)<<make_pair(1,1)
	);
	while(!fit_mc.ConcentratedInOnePoint()){
		fit_mc.Iterate();
		printf("%i iterations \r",fit_mc.iteration_count());
	}
	printf("\n done. Chi^2=%f\n",fit_mc.Optimality());
	for(int i=0,n=fit_mc.ParamCount(); i<n;i++)
		printf("par%i=%f\n",i,fit_mc[i]);
	
	{char* olddir=getcwd(NULL,0);
		chdir(outpath.c_str());
		ofstream data;
		data.open("output.mc.data.txt");
		if(data.is_open()){
			for(int i=0;i<points_mc->Count();i++)
				data<<points_mc->X(i)[0]<<" "
				<<points_mc->Y(i)<<" "
				<<points_mc->X_w(i)[0]<<" "
				<<points_mc->W(i)<<"\n";
			data.close();
		}
		ofstream out;
		out.open("output.mc.func.txt");
		if(out.is_open()){
			for(double x=points_mc->X(0)[0],
				upper=points_mc->X(points_mc->Count()-1)[0];
			x<=upper; x+=0.00001){
				out<<x<<" "<<fit_mc(ParamSet(x))<<"\n";
			}
			out.close();
		}
		ofstream script;
		script.open(".plotscript.gp");
		if(script.is_open()){
			script << "plot ";
			script << "\"output.mc.data.txt\" using 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars title \"data\"";
			script << ",\\\n";
			script << "\"output.mc.func.txt\" w l title \"total fit\"";
			script << "\n";
			script << "\npause -1";
			script.close();
		}
		system("gnuplot .plotscript.gp");
		chdir(olddir);
	}
}