// this file is distributed under 
// MIT license
#ifndef xoyoyptv
#define xoyoyptv
#include <functional>
#include <vector>
#include <map>
#include <memory>
#include <string>
#include <TH1F.h>
#include <TH2F.h>
#include <WTrack.hh>
namespace TrackAnalyse{
	using namespace std;
	class Binner{
	public:
		Binner(double f, double t,unsigned int b);
		Binner(const Binner&source);
		Binner&operator=(const Binner&source);
		virtual ~Binner();
		double left()const;
		double right()const;
		size_t count()const;
		double bin_width() const;
		double bin_center(size_t i) const;
		bool FindBinIndex(size_t&output,double x)const;
	private:
		double from;
		double to;
		size_t bins;
		void CheckCorrectness()const;
	};
	class AHist1D{
	public:
		AHist1D(const string&dir,const string&name,const Binner&x);
		virtual ~AHist1D();
		void Accept(double x)const;
	private:
		TH1F *hist;
	};
	class ASetOfHists1D{
	public:
		ASetOfHists1D(const string&dir,const string&name,const Binner&z,const Binner&x);
		virtual ~ASetOfHists1D();
		void Accept(double z, double x)const;
	private:
		Binner Z;
		TH1F *All,*OutOfBorder;
		vector<TH1F*> Bins;
	};
	class AHist2D{
	public:
		AHist2D(const string&dir,const string&name,const Binner&x,const Binner&y);
		virtual ~AHist2D();
		void Accept(double x,double y)const;
	private:
		TH2F *hist;
	};
	class ASetOfHists2D{
	public:
		ASetOfHists2D(const string&dir,const string&name,const Binner&binning,const Binner&x,const Binner&y);
		virtual ~ASetOfHists2D();
		void Accept(double z,double x, double y)const;
	private:
		Binner Z;
		TH2F *All,*OutOfBorder;
		vector<TH2F*> Bins;
	};


	template<typename... Args>
	class Axis:public Binner{
	private:
		function<double(Args...)> func;
	public:
		Axis():Binner(0,1,1){func=[](Args...){};}
		Axis(function<double(Args...)> v,const Binner&B)
			:Binner(B){func=v;}
		Axis(function<double(Args...)> v,double f, double t,size_t b)
			:Binner(f,t,b){func=v;}
		Axis(const Axis&source):Binner(source){func=source.func;}
		Axis&operator=(const Axis&source){
			Binner::operator=(source);
			func=source.func;
			return *this;
		}
		~Axis(){}
		function<double(Args...)> FUNC()const{return func;}
		double operator()(Args... args)const{return func(args...);}
		bool FindBinNumber(size_t&out,Args... args)const{return FindBinIndex(out,func(args...));}
	};

	template<typename... Args>
	class IProcess{
	public:
		virtual bool Process(Args...)const=0;
		virtual ~IProcess(){}
	};

	template<typename... Args>
	class Hist1D:public IProcess<Args...>,protected AHist1D{
	private:
		function<double(Args...)> X;
	public:
		Hist1D(string&&dir,string&&name,const Axis<Args...>&x):AHist1D(dir,name,x){X=x.FUNC();}
		virtual ~Hist1D(){}
		virtual bool Process(Args...args)const override{
			Accept(X(args...));
			return true;
		};
	};

	template<typename... Args>
	class Hist2D:public IProcess<Args...>,protected AHist2D{
	private:
		function<double(Args...)> X,Y;
	public:
		Hist2D(string&&dir,string&&name,const Axis<Args...>&x,const Axis<Args...>&y):AHist2D(dir,name,x,y){X=x.FUNC();Y=y.FUNC();}
		virtual ~Hist2D(){}
		virtual bool Process(Args...args)const override{
			Accept(X(args...),Y(args...));
			return true;
		};
	};

	template<typename... Args>
	class SetOfHists1D:public IProcess<Args...>,protected ASetOfHists1D{
	private:
		function<double(Args...)> Z,X;
	public:
		SetOfHists1D(string&&dir,string&&name,const Axis<Args...>&z,const Axis<Args...>&x):ASetOfHists1D(dir,name,z,x){Z=z.FUNC();X=x.FUNC();}
		virtual ~SetOfHists1D(){}
		virtual bool Process(Args...args)const override{
			Accept(Z(args...),X(args...));
			return true;
		};
	};

	template<typename... Args>
	class SetOfHists2D:public IProcess<Args...>,protected ASetOfHists2D{
	private:
		function<double(Args...)> Z,X,Y;
	public:
		SetOfHists2D(string&&dir,string&&name,const Axis<Args...>&z,const Axis<Args...>&x,const Axis<Args...>&y):ASetOfHists2D(dir,name,z,x,y){Z=z.FUNC();X=x.FUNC();Y=y.FUNC();}
		virtual ~SetOfHists2D(){}
		virtual bool Process(Args...args)const override{
			Accept(Z(args...),X(args...),Y(args...));
			return true;
		};
	};

	template<typename... Args>
	class Function:public IProcess<Args...>{
	private:
		function<bool(Args...)> func;
	public:
		Function(function<bool(Args...)> f){func=f;}
		virtual ~Function(){}
		virtual bool Process(Args...args)const override{
			return func(args...);
		}
	};

	template<typename... Args>
	class AbstractSet{
	private:
		vector<shared_ptr<IProcess<Args...>>> m_chain;
	protected:
		AbstractSet(){}
	public:
		virtual ~AbstractSet(){}
		void Add(shared_ptr<IProcess<Args...>>element){
			m_chain.push_back(element);
		}
		void Add(function<bool(Args...)>f){
			Add(make_shared<Function<Args...>>(f));
		}
	protected:
		vector<shared_ptr<IProcess<Args...>>>&allset()const{
			return const_cast<vector<shared_ptr<IProcess<Args...>>>&>(m_chain);
		};
	};
	template<class Container,typename... Args>
	class ProcessAll:public Container{
	public:
		ProcessAll(){}
		virtual ~ProcessAll(){}
	protected:
		bool Cycle(Args...args)const{
			for(auto p:AbstractSet<Args...>::allset())
				p->Process(args...);
			return true;
		}
	};
	
	template<class Container,typename... Args>
	class Chain:public Container{
	public:
		Chain(){}
		virtual ~Chain(){}
	protected:
		bool Cycle(Args...args)const{
			for(auto p:AbstractSet<Args...>::allset())
				if(!p->Process(args...))return false;
			return true;
		}
	};
	
	template<class Container,typename... Args>
	class LAnd:public Container{
	public:
		LAnd(){}
		virtual ~LAnd(){}
	protected:
		bool Cycle(Args...args)const{
			bool res=true;
			for(auto p:AbstractSet<Args...>::allset())
				res&=p->Process(args...);
			return res;
		}
	};
	
	template<class Container,typename... Args>
	class LOr:public Container{
	public:
		LOr(){}
		virtual ~LOr(){}
	protected:
		bool Cycle(Args...args)const{
			bool res=false;
			for(auto p:AbstractSet<Args...>::allset())
				res|=p->Process(args...);
			return res;
		}
	};
	
	//WTrack cannot be transfered as const
	//Because all it's methods are not const
	//I have to deal with it
	typedef IProcess<> ParamlessProcess;
	typedef IProcess<WTrack&> ITrackProcess;
	class TrackProcesses:public ITrackProcess,public ProcessAll<AbstractSet<WTrack&>,WTrack&>{
	public:
		virtual bool Process(WTrack&T)const override{return Cycle(T);}
	};
	template<typename element>
	inline shared_ptr<TrackProcesses> operator<<(shared_ptr<TrackProcesses> S,element e){
		dynamic_pointer_cast<AbstractSet<WTrack&>>(S)->Add(e);
		return S;
	}
	typedef vector<double> Results;
	typedef IProcess<const Results&> IResAnalysis;
	class AnalysisSet:public IResAnalysis,public ProcessAll<AbstractSet<const Results&>,const Results&>{
	public:
		virtual bool Process(const Results&P)const override{return Cycle(P);}
	};
	template<typename element>
	inline shared_ptr<AnalysisSet> operator<<(shared_ptr<AnalysisSet> S,element e){
		dynamic_pointer_cast<AbstractSet<const Results&>>(S)->Add(e);
		return S;
	}
	typedef IProcess<WTrack&,Results&> ICalculation;
	class CalculationSet:public AbstractSet<WTrack&,Results&>{
	public:
		void AddC(shared_ptr<ICalculation>element){Add(element);}
		void AddC(function<bool(WTrack&,Results&)>f){Add(f);}
		void AddC(shared_ptr<ParamlessProcess>);
		void AddC(shared_ptr<ITrackProcess>);
		void AddC(shared_ptr<IResAnalysis>);
	};
	template<typename element>
	inline shared_ptr<CalculationSet> operator<<(shared_ptr<CalculationSet> S,element e){
		S->AddC(e);
		return S;
	}
	class TrackCalculations:public ITrackProcess,public ProcessAll<CalculationSet,WTrack&,Results&>{
	public:
		virtual bool Process(WTrack&T)const override;
	};
	class TrackCalcChain:public ITrackProcess,public Chain<CalculationSet,WTrack&,Results&>{
	public:
		virtual bool Process(WTrack&T)const override;
	};
	class Calculations:public ICalculation,public ProcessAll<CalculationSet,WTrack&,Results&>{
	public:
		virtual bool Process(WTrack&T,Results&P)const override{return Cycle(T,P);}
	};
	class CalcChain:public ICalculation,public Chain<CalculationSet,WTrack&,Results&>{
	public:
		virtual bool Process(WTrack&T,Results&P)const override{return Cycle(T,P);}
	};
	class CalcAnd:public ICalculation,public LAnd<CalculationSet,WTrack&,Results&>{
	public:
		virtual bool Process(WTrack&T,Results&P)const override{return Cycle(T,P);}
	};
	class CalcOr:public ICalculation,public LOr<CalculationSet,WTrack&,Results&>{
	public:
		virtual bool Process(WTrack&T,Results&P)const override{return Cycle(T,P);}
	};
}
#endif 