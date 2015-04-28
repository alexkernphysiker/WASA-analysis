#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <memory>
#include <unistd.h>
#include <replace.cc>
#include "gnuplot.h"
using namespace std;
class Plotter{
	vector<string> lines;
	unsigned int counter;
	string outpath;
public:
	Plotter(){
		counter=0;
		outpath="*";
	}
	~Plotter(){
		ofstream script;
		string name="script.gnuplot";
		script.open((outpath+"/"+name).c_str());
		if(script.is_open()){
			for(string line:lines)
				script<<line<<"\n";
			script << "\npause -1";
			script.close();
			name=string("gnuplot ")+name;
			string old=getcwd(NULL,0);
			chdir(outpath.c_str());
			system(name.c_str());
			chdir(old.c_str());
		}
		
	}
	static Plotter &Instance(){
		static Plotter m_instance;
		return m_instance;
	}
	void SetOutput(string out){
		if(outpath!="*")
			throw;
		outpath=out;
	}
	string OutPath(){
		if(outpath=="*")
			throw;
		return outpath;
	}
	string GetTerminal(){
		counter++;
		return string("set terminal pngcairo size 1024,868 enhanced monochrome font 'Verdana,18'\nset output '")+to_string(counter)+".png'";
	}
	Plotter &operator<<(string line){
		lines.push_back(line);
		return *this;
	}
};
void SetPlotOutput(string outpath){
	Plotter::Instance().SetOutput(outpath);
}
Plot::Plot(){
	operator<<(Plotter::Instance().GetTerminal());
}
Plot& Plot::operator<<(string line){
	lines.push_back(line);
}
Plot::~Plot(){
	for(string line:lines)
		Plotter::Instance()<<line;
	for(int i=0,n=plots.size();i<n;i++){
		string line=plots[i];
		if(i==0)
			line="plot "+line;
		if(i<(n-1))
			line+=",\\";
		Plotter::Instance()<<line;
	}
}
#define FILENAME(name) ReplaceAll(ReplaceAll(name," ","_"),"=","_")+".txt"
Plot& Plot::Points(string name, shared_ptr<Genetic::FitPoints> points){
	ofstream data;
	string filename=FILENAME(name);
	data.open((Plotter::Instance().OutPath()+"/"+filename).c_str());
	if(data.is_open()){
		for(auto p:*points)
			data<<p.X[0]<<" "<<p.y<<" "<<p.WX[0]<<" "<<p.wy<<"\n";
		data.close();
		plots.push_back("\""+filename+"\" using 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars title \""+name+"\"");
	}
}
Plot& Plot::Points(string name, LinearInterpolation<double>& points){
	ofstream data;
	string filename=FILENAME(name);
	data.open((Plotter::Instance().OutPath()+"/"+filename).c_str());
	if(data.is_open()){
		for(auto p:points)
			data<<p.first<<" "<<p.second<<"\n";
		data.close();
		plots.push_back("\""+filename+"\" using 1:2 title \""+name+"\"");
	}
}
Plot& Plot::Points(string name, LinearInterpolation<double>& points,function<double(double)> error){
	ofstream data;
	string filename=FILENAME(name);
	data.open((Plotter::Instance().OutPath()+"/"+filename).c_str());
	if(data.is_open()){
		for(auto p:points)
			data<<p.first<<" "<<p.second<<" "<<error(p.first)<<"\n";
		data.close();
		plots.push_back("\""+filename+"\" using 1:2:($2-$3):($2+$3) with yerrorbars title \""+name+"\"");
	}
}
Plot& Plot::Function(string name,function<double(double)> func,double from,double to,double step){
	ofstream out;
	string filename=FILENAME(name);
	out.open((Plotter::Instance().OutPath()+"/"+filename).c_str());
	if(out.is_open()){
		for(double x=from;x<=to;x+=step){
			out<<x<<" "<<func(x)<<"\n";
		}
		out.close();
		plots.push_back("\""+filename+"\" w l title \""+name+"\"");
	}
}
Plot& Plot::Function(string name, function< double(double) > func, function< double(double) > error, double from, double to, double step){
	ofstream out;
	string filename=FILENAME(name);
	out.open((Plotter::Instance().OutPath()+"/"+filename).c_str());
	if(out.is_open()){
		for(double x=from;x<=to;x+=step){
			out<<x<<" "<<func(x)<<" "<<error(x)<<"\n";
		}
		out.close();
		plots.push_back("\""+filename+"\" using 1:2:($2-$3):($2+$3) with yerrorbars title \""+name+"\"");
	}
}
