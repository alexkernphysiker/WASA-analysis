// this file is distributed under 
// GPL v 3.0 license
#include <fstream>
#include "read_simulation.h"
using namespace std;
using namespace Genetic;
shared_ptr<FitPoints> ReadWeightedFrom2D(string name,double from, double to, unsigned int bins){
	auto res=make_shared<FitPoints>();
	FitPoints::Point p[bins][bins];
	double step=(to-from)/double(bins);
	for(unsigned int i=0;i<bins;i++){
		double x=from+step*double(i);
		for(unsigned int j=0;j<bins;j++){
			double y=from+step*double(j);
			p[i][j].X=ParamSet(x);
			p[i][j].y=y;
			p[i][j].wy=0;
		}
	}
	ifstream file;
	file.open(name.c_str());
	if(file.is_open()){
		while(!file.eof()){
			double x,y;
			file>>x>>y;
			int i=int((x-from)/step);
			int j=int((y-from)/step);
			if((i>=0)&&(j>=0)&&(i<bins)&&(j<bins))
				p[i][j].wy+=1.0;
		}
		file.close();
	}
	for(unsigned int i=0;i<bins;i++)
		for(unsigned int j=0;j<bins;j++)
			if(p[i][j].wy>0)
				res<<p[i][j];
	return res;
}
