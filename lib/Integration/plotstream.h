// this file is distributed under 
// MIT license
#ifndef IPZWCLTH
# define IPZWCLTH
#include <string>
#include <utility>
#include <functional>
#include <Genetic/fit.h>
#include <gnuplot_wrap.h>
namespace PlotStream{
	using namespace std;
	using namespace Genetic;
	typedef PlotPoints<double,vector<pair<double,double>>> PlotEngine;
	class AbstractPlotStream{
	protected:
		AbstractPlotStream(string&&name,shared_ptr<PlotEngine>plot);
		AbstractPlotStream(string&&name);
		virtual ~AbstractPlotStream();
	public:
		AbstractPlotStream&operator<<(const ParamSet&P);
		AbstractPlotStream&operator<<(ParamSet&&P);
		AbstractPlotStream&operator<<(const Point&P);
		AbstractPlotStream&operator<<(Point&&P);
		shared_ptr<PlotEngine>Plot();
	protected:
		virtual void ProcessPoint(const ParamSet&P)=0;
		string&&Name()const;
	private:
		string m_name;
		shared_ptr<PlotEngine> m_plot;
	};
	class SimplePlotStream:public AbstractPlotStream{
	public:
		SimplePlotStream(string&& name,pair<size_t,size_t>&&indexes);
		SimplePlotStream(string&& name,pair<size_t,size_t>&&indexes,shared_ptr<PlotEngine> plot);
		virtual ~SimplePlotStream();
		typedef function<double(double)> func;
		SimplePlotStream&AddFunc(func f);
	protected:
		virtual void ProcessPoint(const ParamSet& P)override;
	private:
		pair<size_t,size_t> m_indexes;
		vector<pair<double,double>> m_data;
		pair<double,double> m_xrange;
		vector<func> m_funcs;
	};
	class Binner:public AbstractPlotStream{
	public:
		struct binparam{
			size_t param_index;
			double from,to;
			size_t bin_count;
			binparam(size_t ind,double a,double b,size_t cnt);
			binparam(const binparam&source);
		};
		Binner(string&& name,const binparam& binning);
		Binner(string&& name,binparam&&binning);
		virtual ~Binner();
		typedef function<shared_ptr<SimplePlotStream>(size_t)>Creator1;
		typedef function<shared_ptr<SimplePlotStream>(size_t,string&&)>Creator2;
		typedef function<shared_ptr<SimplePlotStream>(size_t,string&&,shared_ptr<PlotEngine>)>Creator3;
		Binner&Fill(Creator1 func);
		Binner&Fill(Creator2 func);
		Binner&Fill(Creator3 func);
		typedef function<double(double b,double x)> func;
		Binner&AddFunc(func f);
	protected:
		virtual void ProcessPoint(const ParamSet& P)override;
		bool FindBinIndex(const ParamSet& P,size_t&res)const;
		double BinHW()const;
		double BinCenter(size_t i)const;
		size_t BinCount()const;
	private:
		void CheckCorrectness()const;
		binparam m_binning;
		vector<shared_ptr<SimplePlotStream>> m_streams;
	};
};
#endif