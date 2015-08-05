// this file is distributed under 
// GPL v 3.0 license
#include <fstream>
#include "read_simulation.h"
using namespace std;
using namespace Genetic;
shared_ptr< FitPoints > ReadWeightedFrom2D(
	string name, 
	double fromX, double toX, unsigned int binsX, 
	double fromY, double toY, unsigned int binsY, 
	function< bool > processing
){
	auto res=make_shared<FitPoints>();
	FitPoints::Point p[binsX][binsY];
	double stepX=(toX-fromX)/double(binsX);
	double stepY=(toY-fromY)/double(binsY);
	for(unsigned int i=0;i<binsX;i++){
		double x=fromX+stepX*double(i);
		for(unsigned int j=0;j<binsY;j++){
			double y=fromY+stepY*double(j);
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
			if(processing(x,y)){
				int i=int((x-fromX)/stepX);
				int j=int((y-fromY)/stepY);
				if((i>=0)&&(j>=0)&&(i<binsX)&&(j<binsY))
					p[i][j].wy+=1.0;
			}
		}
		file.close();
	}
	for(unsigned int i=0;i<binsX;i++)
		for(unsigned int j=0;j<binsY;j++)
			if(p[i][j].wy>0)
				res<<p[i][j];
	return res;
}
