// this file is distributed under 
// MIT license
#ifndef IPZWCLTH
# define IPZWCLTH
#include <string>
#include <utility>
#include <Genetic/fit.h>
namespace PlotStream{
	using namespace std;
	using namespace Genetic;
	class AbstractPlotStream{
	protected:
		AbstractPlotStream(string&&name);
	public:
		AbstractPlotStream&operator<<(const ParamSet&P);
		AbstractPlotStream&operator<<(ParamSet&&P);
		AbstractPlotStream&operator<<(const Point&P);
		AbstractPlotStream&operator<<(Point&&P);
	protected:
		virtual void ProcessPoint(const ParamSet&P)=0;
		string&Name()const;
	private:
		string m_name;
	};
	class SimpleStream:public AbstractPlotStream{
	public:
		SimpleStream(string&& name,pair<size_t,size_t>&&indexes);
		virtual ~SimpleStream();
	protected:
		virtual void ProcessPoint(const ParamSet& P)override;
	private:
		vector<pair<double,double>> m_data;
	};
};
#endif