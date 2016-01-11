// this file is distributed under 
// MIT license
#include <math_h/exception_math_h.h>
#include "plotstream.h"
namespace PlotStream {
	using namespace std;
	using namespace Genetic;
	ParamSet GM(const Point&p){
		ParamSet res=p.X();
		res<<p.y();
		return res;
	}
AbstractPlotStream::AbstractPlotStream(string&& name, shared_ptr<PlotEngine> plot){
	m_plot=plot;
	m_name=name;
}
AbstractPlotStream::AbstractPlotStream(string&& name):
	AbstractPlotStream(static_cast<string&&>(name),make_shared<PlotEngine>()){}
string&&AbstractPlotStream::Name() const{return const_cast<string&&>(m_name);}
shared_ptr< PlotEngine > AbstractPlotStream::Plot(){return m_plot;}
AbstractPlotStream::~AbstractPlotStream(){}
AbstractPlotStream& AbstractPlotStream::operator<<(const ParamSet& P){
	ProcessPoint(P);
	return *this;
}
AbstractPlotStream& AbstractPlotStream::operator<<(ParamSet&& P){return operator<<(P);}
AbstractPlotStream& AbstractPlotStream::operator<<(const Point& P){return operator<<(GM(P));}
AbstractPlotStream& AbstractPlotStream::operator<<(Point&& P){return operator<<(GM(P));}
SimplePlotStream::SimplePlotStream(string&&name, pair<size_t,size_t>&&indexes,shared_ptr<PlotEngine>plot)
	:AbstractPlotStream(static_cast<string&&>(name), plot){m_indexes=indexes;}

SimplePlotStream::SimplePlotStream(std::string&& name,pair<size_t,size_t>&& indexes)
	:AbstractPlotStream(static_cast<string&&>(name)){m_indexes=indexes;}
SimplePlotStream::~SimplePlotStream(){
	Plot()->WithoutErrors(Name(),m_data);
}
void SimplePlotStream::ProcessPoint(const ParamSet& P){
	m_data.push_back(make_pair(P[m_indexes.first],P[m_indexes.first]));
}

Binner::binparam::binparam(size_t ind, double a, double b, size_t cnt)
	:param_index(ind),from(a),to(b),bin_count(cnt){}
Binner::binparam::binparam(const Binner::binparam& source)
	:param_index(source.param_index)
	,from(source.from),to(source.to)
	,bin_count(source.bin_count){}

Binner::Binner(string&&name,const Binner::binparam&binning)
:AbstractPlotStream(static_cast<string&&>(name)),m_binning(binning){
	CheckCorrectness();
}
Binner::Binner(string&& name, Binner::binparam&& binning)
:AbstractPlotStream(static_cast<string&&>(name)),m_binning(binning){
	CheckCorrectness();
}
void Binner::CheckCorrectness()const{
	if(m_binning.to<m_binning.from)
		throw math_h_error<decltype(*this)>("wrong binning ranges");
	if(0==m_binning.bin_count)
		throw math_h_error<decltype(*this)>("there cannot be zero bins");
}
Binner::~Binner(){}
size_t Binner::BinCount()const{return m_binning.bin_count;}
double Binner::BinHW() const{
	return (m_binning.to-m_binning.from)/double(m_binning.bin_count*2);
}
double Binner::BinCenter(size_t i) const{
	if(i>=m_binning.bin_count)
		throw math_h_error<decltype(*this)>("bin range check error");
	return m_binning.from+BinHW()*double(1+2*i);
}
bool Binner::FindBinIndex(const ParamSet& P, size_t& res) const{
	double x=P[m_binning.param_index],delta=BinHW();
	for(size_t i=0,n=BinCount();i<n;i++)
		if(
			(x>=(BinCenter(i)-delta))&&
			(x<(BinCenter(i)+delta))
		){
			res=i;
			return true;
		}
	return false;
}
void Binner::Fill(Binner::Creator1 func){
	m_streams.clear();
	for(size_t i=0,n=BinCount();i<n;i++)
		m_streams.push_back(func(i));
}
void Binner::Fill(Binner::Creator2 func){
	m_streams.clear();
	for(size_t i=0,n=BinCount();i<n;i++)
		m_streams.push_back(func(i,Name()));
}
void Binner::Fill(Binner::Creator3 func){
	m_streams.clear();
	for(size_t i=0,n=BinCount();i<n;i++)
		m_streams.push_back(func(i,Name(),Plot()));
}
void Binner::ProcessPoint(const ParamSet& P){
	if(m_streams.size()==0)
		throw math_h_error<decltype(*this)>("Attempt to use not filled binner");
	size_t index=0;
	if(FindBinIndex(P,index))
		m_streams[index]->operator<<(P);
}
}