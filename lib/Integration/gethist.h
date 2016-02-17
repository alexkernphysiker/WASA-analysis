// this file is distributed under 
// MIT license
#ifndef DAIHADDG
#define DAIHADDG
#include <string>
#include <functional>
#include <vector>
#include <math_h/hist.h>
namespace ROOT_data{
	using namespace std;
	using namespace MathTemplates;
	double PresentRunsAmountRatio(string&&reaction);
	enum histsource{MC,DATA};
	hist<double> Hist(histsource src, const string&reaction, const vector<string>&path,const string&histname);
	hist<double> Hist(histsource src,const string&reaction,const vector<string>&path,string&&histname);
	hist<double> Hist(histsource src,string&&reaction,const vector<string>&path,string&&histname);
	hist<double> Hist(histsource src,string&&reaction,vector<string>&&path,string&&histname);
};
#endif