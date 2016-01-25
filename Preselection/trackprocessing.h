// this file is distributed under 
// MIT license
#ifndef xoyoyptv
#define xoyoyptv
#include <functional>
#include <vector>
#include <memory>
#include <string>
#include "analysis.h"
namespace TrackAnalyse{
	using namespace std;
	//WTrack cannot be transfered as const because
	//it does not contain const methods 
	//(even those ones that really should be)
	typedef function<double()> ValueIndependent;
	typedef function<double(WTrack&)> ValueTrackDependent;
	typedef function<double(const vector<double>&)> ValueParamDependent;
	typedef function<double(WTrack&,const vector<double>&)> ValueTrackParamDependent;
	typedef function<bool()> ConditionIndependent;
	typedef function<bool(WTrack&)> ConditionTrackDependent;
	typedef function<bool(const vector<double>&)> ConditionParamDependent;
	typedef function<bool(WTrack&,const vector<double>&)> ConditionTrackParamDependent;
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
	class ITrackParamProcess{
	public:
		virtual bool Process(WTrack&,vector<double>&)const=0;
		virtual ~ITrackParamProcess(){}
	};
	class ITrackParamAnalyse:public ITrackParamProcess{
	public:
		virtual bool Process(WTrack&,vector<double>&)const final;
		virtual ~ITrackParamAnalyse(){}
	protected:
		virtual void Analyse(WTrack&,const vector<double>&)const =0;
	};
	
	class Hist1D:public ITrackParamAnalyse{
	public:
		Hist1D(string&&name,const Axis&x);
		virtual ~Hist1D();
	protected:
		virtual void Analyse(WTrack&,const vector<double>&)const override;
	private:
		TH1F *hist;
		Axis X;
	};
	class SetOfHists1D:public ITrackParamAnalyse{
	public:
		SetOfHists1D(string&&name,const Axis&binning,const Axis&x);
		virtual ~SetOfHists1D();
	protected:
		virtual void Analyse(WTrack&,const vector<double>&)const override;
	private:
		Axis Z,X;
		TH1F *All,*OutOfBorder;
		vector<TH1F*> Bins;
	};

	class Hist2D:public ITrackParamAnalyse{
	public:
		Hist2D(string&&name,const Axis&x,const Axis&y);
		virtual ~Hist2D();
	protected:
		virtual void Analyse(WTrack&,const vector<double>&)const override;
	private:
		TH2F *hist;
		Axis X,Y;
	};
	class SetOfHists2D:public ITrackParamAnalyse{
	public:
		SetOfHists2D(string&&name,const Axis&binning,const Axis&x,const Axis&y);
		virtual ~SetOfHists2D();
	protected:
		virtual void Analyse(WTrack&,const vector<double>&)const override;
	private:
		Axis Z,X,Y;
		TH2F *All,*OutOfBorder;
		vector<TH2F*> Bins;
	};
	
	class Condition:public ITrackParamProcess{
	public:
		Condition(ConditionTrackParamDependent func);
		Condition(ConditionTrackDependent func);
		Condition(ConditionParamDependent func);
		Condition(ConditionIndependent func);
		virtual ~Condition();
		virtual bool Process(WTrack&,vector<double>&)const override;
	private:
		ConditionTrackParamDependent condition;
	};
	class Parameter:public ITrackParamProcess{
	public:
		Parameter(ValueTrackParamDependent f);
		Parameter(ValueParamDependent f);
		Parameter(ValueTrackDependent f);
		Parameter(ValueIndependent f);
		virtual ~Parameter();
		virtual bool Process(WTrack&,vector<double>&)const override;
	private:
		ValueTrackParamDependent func;
	};
	class AbstractChain:public ITrackParamProcess{
	protected:
		AbstractChain();
	public:
		virtual ~AbstractChain();
		AbstractChain&operator<<(shared_ptr<ITrackParamProcess>element);
	protected:
		void Cycle(function<void(bool)>,WTrack&,vector<double>&)const;
	private:
		vector<shared_ptr<ITrackParamProcess>> m_chain;
	};
	inline shared_ptr<AbstractChain>operator<<(shared_ptr<AbstractChain>ch,shared_ptr<ITrackParamProcess>v){
		ch->operator<<(v);
		return ch;
	}
	inline shared_ptr<AbstractChain>operator<<(shared_ptr<AbstractChain>ch,ValueTrackParamDependent f){
		return ch<<make_shared<Parameter>(f);
	}
	inline shared_ptr<AbstractChain>operator<<(shared_ptr<AbstractChain>ch,ValueParamDependent f){
		return ch<<make_shared<Parameter>(f);
	}
	inline shared_ptr<AbstractChain>operator<<(shared_ptr<AbstractChain>ch,ValueTrackDependent f){
		return ch<<make_shared<Parameter>(f);
	}
	inline shared_ptr<AbstractChain>operator<<(shared_ptr<AbstractChain>ch,ValueIndependent f){
		return ch<<make_shared<Parameter>(f);
	}
	inline shared_ptr<AbstractChain>operator<<(shared_ptr<AbstractChain>ch,ConditionTrackParamDependent f){
		return ch<<make_shared<Condition>(f);
	}
	inline shared_ptr<AbstractChain>operator<<(shared_ptr<AbstractChain>ch,ConditionTrackDependent f){
		return ch<<make_shared<Condition>(f);
	}
	inline shared_ptr<AbstractChain>operator<<(shared_ptr<AbstractChain>ch,ConditionParamDependent f){
		return ch<<make_shared<Condition>(f);
	}
	inline shared_ptr<AbstractChain>operator<<(shared_ptr<AbstractChain>ch,ConditionIndependent f){
		return ch<<make_shared<Condition>(f);
	}
	class ChainAnd:public AbstractChain{
	public:
		ChainAnd(){}
		virtual ~ChainAnd(){}
		virtual bool Process(WTrack&,vector<double>&)const override;
	};
	class ChainOr:public AbstractChain{
	public:
		ChainOr(){}
		virtual ~ChainOr(){}
		virtual bool Process(WTrack&,vector<double>&)const override;
	};
	class TrackProcess{
	public:
		TrackProcess(){}
		~TrackProcess(){}
		TrackProcess&operator<<(shared_ptr<ITrackParamProcess>element);
		void Process(WTrack&T)const;
	private:
		vector<shared_ptr<ITrackParamProcess>> m_proc;
	};
}
#endif 