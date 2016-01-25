// this file is distributed under 
// MIT license
#ifndef xoyoyptv
#define xoyoyptv
#include <functional>
#include <vector>
#include <string>
#include "analysis.h"
//WTrack cannot be transfered as const because
//it does not contain const methods 
//(even those ones that really should be)
typedef std::function<double()> ValueIndependent;
typedef std::function<double(WTrack&)> ValueTrackDependent;
typedef std::function<double(const std::vector<double>&)> ValueParamDependent;
typedef std::function<double(WTrack&,const std::vector<double>&)> ValueTrackParamDependent;
typedef std::function<bool()> ConditionIndependent;
typedef std::function<bool(WTrack&)> ConditionTrackDependent;
typedef std::function<bool(const std::vector<double>&)> ConditionParamDependent;
typedef std::function<bool(WTrack&,const std::vector<double>&)> ConditionTrackParamDependent;
class Analyser2D;
class TrackConditionSet{
	friend class Analyser2D;
public:
	typedef std::function<bool(WTrack&,std::vector<double>&)> InternalCondition;
	TrackConditionSet(std::string&&name,ValueIndependent distr,int bins,double from,double to);
	TrackConditionSet(std::string&&name,const TrackConditionSet&master);
	virtual ~TrackConditionSet();
	TrackConditionSet&AddParameter(std::string&&name,ValueTrackDependent parameter);
	TrackConditionSet&AddConditions(string&& name,std::vector<TrackConditionSet>&set);
	TrackConditionSet&AddCondition(std::string&&name,ConditionTrackParamDependent condition);
	TrackConditionSet&AddCondition(std::string&&name,ConditionTrackDependent condition);
	TrackConditionSet&AddCondition(std::string&&name,ConditionParamDependent condition);
	TrackConditionSet&AddCondition(std::string&&name,ConditionIndependent condition);
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
		TrackCondition(const std::string&n,InternalCondition delegate,TrackConditionSet*master);
        virtual ~TrackCondition();
		bool Check(WTrack&,std::vector<double>&,double magnitude);
	private:
		InternalCondition condition;
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
	void AcceptEvent(const std::vector<double>&);
private:
	std::string m_name;
	ValueParamDependent m_B;
	TH1F *ForAllA;
	std::vector<TH1F*> A_bin;
	TH1F* reference;
	ValueIndependent m_distr;
};
class Axis{
public:
	Axis(ValueTrackParamDependent v,double f, double t,unsigned int b);
	Axis(ValueTrackDependent v,double f, double t,unsigned int b);
	Axis(ValueParamDependent v,double f, double t,unsigned int b);
	Axis(ValueIndependent v,double f, double t,unsigned int b);
	Axis(const Axis&source);
	~Axis();
	double left()const;
	double right()const;
	unsigned int count()const;
	double getvalue(WTrack&T,const std::vector<double>&P)const;
	ValueTrackParamDependent valuegetter()const;
	double bin_width()const;
	double bin_center(size_t i)const;
	bool FindBinIndex(unsigned int&output,WTrack&T,const std::vector<double>&P)const;
private:
	void CheckCorrectness()const;
	ValueTrackParamDependent value;
	double from;
	double to;
	unsigned int bins;
};
class Debug2DSpectraSet{
public:
	Debug2DSpectraSet(string&&name);
	virtual ~Debug2DSpectraSet();
	void CatchState(WTrack&track,const vector<double>&P);
	typedef std::pair<double,double> point;
	typedef std::function<point(WTrack&,const vector<double>&)> Process;
	typedef std::pair<TH2F*,Process> Item;
	void Add(string&&name,const Axis&X,const Axis&Y);
	void Add(string&&name,Axis&&X,Axis&&Y);
private:
	std::string m_name;
	std::vector<Item> jobs;
	Process Create(ValueTrackParamDependent x,ValueTrackParamDependent y);
};
class Debug2DSpectraBined{
public:
	Debug2DSpectraBined(string&&name,const Axis&x,const Axis&y,const Axis&z);
	Debug2DSpectraBined(string&&name,Axis&&x,Axis&&y,Axis&&z);
	virtual ~Debug2DSpectraBined();
	void CatchState(WTrack&track,const vector<double>&P);
private:
	std::vector<TH2F*> sp2;
	Axis X,Y,Z;
};
#endif 