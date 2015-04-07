#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <memory>
#include <unistd.h>
#include "gnuplot.h"
using namespace std;
Plot::Plot(string out,string script){
	outpath=out;
	scriptname=script;
}
Plot::~Plot(){
	ofstream script;
	script.open((outpath+"/"+scriptname+".gnuplot").c_str());
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
}
Plot& Plot::Points(string file, shared_ptr< Fit::FitPoints> points){
	ofstream data;
	data.open((outpath+"/"+file+".txt").c_str());
	if(data.is_open()){
		for(auto p:*points)
			data<<p.X[0]<<" "<<p.y<<" "<<p.WX[0]<<" "<<p.wy<<"\n";
		data.close();
		lines.push_back("\""+file+".txt\" using 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars title \""+file+"\"");
	}
}

Plot& Plot::Points(string file, hist& points){
	ofstream data;
	data.open((outpath+"/"+file+".txt").c_str());
	if(data.is_open()){
		for(hist::point p:points)
			data<<p.x<<" "<<p.y<<" "<<p.dx<<" "<<p.dy<<"\n";
		data.close();
		lines.push_back("\""+file+".txt\" using 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars title \""+file+"\"");
	}
}
Plot& Plot::Points(string file, LinearInterpolation<double>& points){
	ofstream data;
	data.open((outpath+"/"+file+".txt").c_str());
	if(data.is_open()){
		for(auto p:points)
			data<<p.first<<" "<<p.second<<"\n";
		data.close();
		lines.push_back("\""+file+".txt\" using 1:2 title \""+file+"\"");
	}
}
Plot& Plot::Points(string file, LinearInterpolation<double>& points,function<double(double)> error){
	ofstream data;
	data.open((outpath+"/"+file+".txt").c_str());
	if(data.is_open()){
		for(auto p:points)
			data<<p.first<<" "<<p.second<<" "<<error(p.first)<<"\n";
		data.close();
		lines.push_back("\""+file+".txt\" using 1:2:($2-$3):($2+$3) with yerrorbars title \""+file+"\"");
	}
}
Plot& Plot::Function(string file,function<double(double)> func,double from,double to,double step){
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
Plot& Plot::Function(string file, function< double(double) > func, function< double(double) > error, double from, double to, double step){
	ofstream out;
	out.open((outpath+"/"+file+".txt").c_str());
	if(out.is_open()){
		for(double x=from;x<=to;x+=step){
			out<<x<<" "<<func(x)<<" "<<error(x)<<"\n";
		}
		out.close();
		lines.push_back("\""+file+".txt\" using 1:2:($2-$3):($2+$3) with yerrorbars title \""+file+"\"");
	}
}
