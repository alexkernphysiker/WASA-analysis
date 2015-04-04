#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <memory>
#include <unistd.h>
#include "gnuplot.h"
using namespace std;
Plot::Plot(std::string out){
	outpath=out;
}
Plot::~Plot(){}
Plot& Plot::Points(std::string file, std::shared_ptr< Fit::FitPointsAbstract > points){
	ofstream data;
	data.open((outpath+"/"+file+".txt").c_str());
	if(data.is_open()){
		for(int i=0;i<points->Count();i++)
			data<<points->X(i)[0]<<" "
			<<points->Y(i)<<" "
			<<points->X_w(i)[0]<<" "
			<<points->W(i)<<"\n";
		data.close();
		lines.push_back("\""+file+".txt\" using 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars title \""+file+"\"");
	}
}

Plot& Plot::Hist(std::string file, hist& points){
	ofstream data;
	data.open((outpath+"/"+file+".txt").c_str());
	if(data.is_open()){
		for(point p:points)
			data<<p.x<<" "<<p.y<<" "<<p.dx<<" "<<p.dy<<"\n";
		data.close();
		lines.push_back("\""+file+".txt\" using 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars title \""+file+"\"");
	}
}
Plot& Plot::Function(std::string file,std::function<double(double)> func,double from,double to,double step){
	ofstream out;
	out.open((outpath+"/"+file+".txt").c_str());
	if(out.is_open()){
		for(double x=from;x<=to;x+=step){
			out<<x<<" "<<func(x)<<"\n";
		}
		out.close();
		lines.push_back("\""+file+".txt\" w l title \""+file+"\"");
	}
}
void Plot::Out(string scriptname,bool show){
	ofstream script;
	script.open((outpath+"/"+scriptname).c_str());
	if(script.is_open()){
		script << "plot ";
		for(int i=0,n=lines.size();i<n;i++){
			script<<lines[i];
			if(i<(n-1))
				script<<",\\";
			script<<"\n";
		}
		script << "\npause -1";
		script.close();
	}
	if(show){
		char* olddir=getcwd(NULL,0);
		chdir(outpath.c_str());
		system((string("gnuplot ")+scriptname).c_str());
		chdir(olddir);
	}
}




