// this file is distributed under 
// MIT license
#include <fstream>
#include <str_get.h>
#include "read_simulation.h"
using namespace std;
using namespace Genetic;
string SimulationDataPath(){
	static string str="";
	if(str=="")str=ENV(PRESEL_DATA)+string("/../Reconstruction/");
	return str;
}
shared_ptr<FitPoints> ReadWeightedFrom2D(
	string name, 
	double fromX, double toX, unsigned int binsX, 
	double fromY, double toY, unsigned int binsY, 
	function<bool(double&,double&)> processing
){
	auto res=make_shared<FitPoints>();
	typedef FitPoints::Point point;
	typedef point* ppoint;
	ppoint *p=new ppoint[binsX];
	for(unsigned int x=0;x<binsX;x++)
		p[x]=new point[binsY];
	double stepX=(toX-fromX)/double(binsX);
	double stepY=(toY-fromY)/double(binsY);
	for(unsigned int i=0;i<binsX;i++){
		double x=fromX+stepX*(double(i)+0.5);
		for(unsigned int j=0;j<binsY;j++){
			double y=fromY+stepY*(double(j)+0.5);
			p[i][j].X<<x;
			p[i][j].y=y;
			p[i][j].wy=0;
		}
	}
	ifstream file;ofstream out;
	file.open(name.c_str());
	out.open((name+".cut.txt").c_str());
	if(file.is_open()&&out.is_open()){
		while(!file.eof()){
			double x,y;
			file>>x>>y;
			if(processing(x,y)){
				int i=int((x-fromX)/stepX);
				int j=int((y-fromY)/stepY);
				if((i>=0)&&(j>=0)&&(i<binsX)&&(j<binsY))
					p[i][j].wy+=1.0;
				out<<x<<" "<<y<<endl;
			}
		}
		file.close();
	}else
		throw math_h_error<ifstream>("Not found");
	for(unsigned int i=0;i<binsX;i++)
		for(unsigned int j=0;j<binsY;j++)
			if(p[i][j].wy>0)
				res<<p[i][j];
	for(unsigned int x=0;x<binsX;x++)
		delete[] p[x];
	delete[] p;
	return res;
}
