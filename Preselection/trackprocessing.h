// this file is distributed under 
// GPL v 3.0 license
#ifndef xoyoyptv
#define xoyoyptv
#include <functional>
#include <vector>
#include <string>
#include "analysis.h"
typedef std::function<double()> ValueIndependent;
typedef std::function<double(WTrack&)> ValueTrackDependent;
typedef std::function<double(std::vector<double>&)> ValueParamDependent;
typedef std::function<bool()> ConditionIndependent;
typedef std::function<bool(WTrack&)> ConditionTrackDependent;
typedef std::function<bool(std::vector<double>&)> ConditionParamDependent;
class Analyser2D;
class TrackConditionSet{
	friend class Analyser2D;
public:
	typedef std::function<bool(WTrack&,std::vector<double>&)> Condition;
	TrackConditionSet(std::string&&name,ValueIndependent distr,int bins,double from,double to);
	TrackConditionSet(std::string&&name,const TrackConditionSet&master);
	virtual ~TrackConditionSet();
	TrackConditionSet&AddParameter(std::string&&name,ValueTrackDependent parameter);
	TrackConditionSet&AddCondition(std::string&&name,Condition condition);
	void ReferenceEvent();
	bool Check(WTrack&track,std::vector<double>&parameters);
protected:
	TH1F* reference;
	ValueIndependent m_distr;
private:
	class TrackCalc{
	public:
		TrackCalc(const std::string&n);
		virtual ~TrackCalc();
		std::string&&Name();
	private:
		std::string name;
	};
	class ParamCalc:public TrackCalc{
	public:
		ParamCalc(const std::string&n,ValueTrackDependent);
		virtual ~ParamCalc();
		double Get(WTrack&);
	private:
		ValueTrackDependent m_delegate;
	};
	class TrackCondition:public TrackCalc{
	public:
		TrackCondition(const std::string&n,Condition delegate,TrackConditionSet*master);
        virtual ~TrackCondition();
		bool Check(WTrack&,std::vector<double>&,double magnitude);
	private:
		Condition condition;
		TH1F* output;
	};
	std::vector<std::shared_ptr<TrackCalc>> calc_procs;
	std::string m_name;
	double m_from,m_to; int m_bins;
	TH1F *beforecut;
};
class Analyser2D{
public:
	Analyser2D(std::string&&name,const TrackConditionSet&);
	void Setup(ValueParamDependent B,int binsB,double fromB,double toB);
	virtual ~Analyser2D();
	void AcceptEvent(std::vector<double>&);
private:
	std::string m_name;
	ValueParamDependent m_B;
	TH1F *ForAllA;
	std::vector<TH1F*> A_bin;
	TH1F* reference;
	ValueIndependent m_distr;
};
struct Axis{
	ValueTrackDependent value;
	double from;
	double to;
	int bins;
};
class Debug2DSpectraSet{
public:
	Debug2DSpectraSet(string&&name);
	virtual ~Debug2DSpectraSet();
	void CatchState(WTrack&track);
	typedef std::pair<double,double> point;
	typedef std::function<point(WTrack&)> Process;
	typedef std::pair<TH2F*,Process> Item;
	void Add(string&&name,const Axis&X,const Axis&Y);
private:
	std::string m_name;
	std::vector<Item> jobs;
	Process Create(ValueTrackDependent x,ValueTrackDependent y);
};
#endif 