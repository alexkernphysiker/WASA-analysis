// this file is distributed under 
// MIT license
#include "TH1F.h"
#include "TH2F.h"
#include <TLorentzVector.h>
#include <CAnalysisModule.hh>
#include <CDataManager.hh>
#include <FDEdep2Ekin.hh>
#include <CCardWDET.hh>
#include <Wasa.hh>
#include <CAnalysisModule.hh>
#include <REventWmcHeader.hh>
#include <REventHeader.hh>
#include <WTrackBank.hh>
#include <WVertexBank.hh>
#include <FDFTHTracks.hh>
#include <CDTracksSimple.hh>
#include "math_h/error.h"
#include "trackprocessing.h"
#include "../config.h"
namespace TrackAnalyse{
	using namespace std;
	using namespace MathTemplates;
	Axis::Axis(ValueTrackParamDependent v, double f, double t, unsigned int b){
		value=v;from=f;to=t;bins=b;
		CheckCorrectness();
	}
	Axis::Axis(ValueTrackDependent v, double f, double t, unsigned int b)
	:Axis([v](WTrack&track,const vector<double>&){return v(track);},f,t,b){}
	Axis::Axis(ValueParamDependent v, double f, double t, unsigned int b)
	:Axis([v](WTrack&,const vector<double>&param){return v(param);},f,t,b){}
	Axis::Axis(ValueIndependent v, double f, double t, unsigned int b)
	:Axis([v](WTrack&,const vector<double>&){return v();},f,t,b){}
	Axis::Axis(const Axis& source):Axis(source.value,source.from,source.to,source.bins){}
	Axis::~Axis(){}
	void Axis::CheckCorrectness() const{
		if(to<=from)
			throw Exception<Axis>("wrong binning ranges");
		if(0==bins)
			throw Exception<Axis>("there cannot be zero bins");
	}
	unsigned int Axis::count() const{return bins;}
	double Axis::left() const{return from;}
	double Axis::right() const{return to;}
	double Axis::bin_width() const{
		if(0==bins)throw Exception<Axis>("count==0");
		return (to-from)/double(bins);
	}
	double Axis::bin_center(size_t i) const{
		if(i>=bins)throw Exception<Axis>("bin range check error");
		return from+bin_width()*(double(i)+0.5);
	}
	double Axis::getvalue(WTrack& T, const vector< double >& P) const{
		return value(T,P);
	}
	ValueTrackParamDependent Axis::valuegetter() const{
		return value;
	}
	bool Axis::FindBinIndex(unsigned int& output, WTrack& T, const vector< double >& P) const{
		double x=getvalue(T,P),delta=bin_width()/2.0;
		if(delta<=0)throw Exception<Axis>("delta<=0");
		for(size_t i=0,n=count();i<n;i++){
			double pos=bin_center(i);
			if((x>=(pos-delta))&&(x<(pos+delta))){
				output=i;
				return true;
			}
		}
		return false;
		
	}
	
	bool ITrackParamAnalyse::Process(WTrack&T, vector<double>&P)const{
		Analyse(T,P);
		return true;
	}
	
	Hist1D::Hist1D(string&& name, const Axis& x):X(x){
		hist=new TH1F(name.c_str(),"",X.count(),X.left(),X.right());
		gHistoManager->Add(hist,"GeneralHistograms");
	}
	Hist1D::~Hist1D(){}
	void Hist1D::Analyse(WTrack&T, const vector<double>&P)const{hist->Fill(X.getvalue(T,P));}
	SetOfHists1D::SetOfHists1D(string&& name, const Axis& binning, const Axis& x):Z(binning),X(x){
		All=new TH1F("ForAllBins","",X.count(),X.left(),X.right());
		gHistoManager->Add(All,name.c_str());
		OutOfBorder=new TH1F("OutOfBins","",X.count(),X.left(),X.right());
		gHistoManager->Add(OutOfBorder,name.c_str());
		for(unsigned int i=0;i<Z.count();i++){
			TH1F*hist=new TH1F(to_string(i).c_str(),"",X.count(),X.left(),X.right());
			gHistoManager->Add(hist,name.c_str());
			Bins.push_back(hist);
		}
	}
	SetOfHists1D::~SetOfHists1D(){}
	void SetOfHists1D::Analyse(WTrack&T, const vector<double>&P)const{
		double x=X.getvalue(T,P);
		unsigned int i=0;
		if(Z.FindBinIndex(i,T,P)){
			All->Fill(x);
			Bins[i]->Fill(x);
		}else{
			OutOfBorder->Fill(x);
		}
	}
	Hist2D::Hist2D(string&& name, const Axis& x, const Axis& y):X(x),Y(y){
		hist=new TH2F(name.c_str(),"",X.count(),X.left(),X.right(),Y.count(),Y.left(),Y.right());
		gHistoManager->Add(hist,"GeneralHistograms");
	}
	Hist2D::~Hist2D(){}
	void Hist2D::Analyse(WTrack&T, const vector<double>&P)const{hist->Fill(X.getvalue(T,P),Y.getvalue(T,P));}
	SetOfHists2D::SetOfHists2D(string&& name, const Axis& binning, const Axis& x, const Axis& y):Z(binning),X(x),Y(y){
		All=new TH2F("ForAllBins","",X.count(),X.left(),X.right(),Y.count(),Y.left(),Y.right());
		gHistoManager->Add(All,name.c_str());
		OutOfBorder=new TH2F("OutOfBins","",X.count(),X.left(),X.right(),Y.count(),Y.left(),Y.right());
		gHistoManager->Add(OutOfBorder,name.c_str());
		for(unsigned int i=0;i<Z.count();i++){
			TH2F*hist=new TH2F(to_string(i).c_str(),"",X.count(),X.left(),X.right(),Y.count(),Y.left(),Y.right());
			gHistoManager->Add(hist,name.c_str());
			Bins.push_back(hist);
		}
	}
	SetOfHists2D::~SetOfHists2D(){}
	void SetOfHists2D::Analyse(WTrack&T, const vector<double>&P)const{
		double x=X.getvalue(T,P);
		double y=Y.getvalue(T,P);
		unsigned int i=0;
		if(Z.FindBinIndex(i,T,P)){
			All->Fill(x,y);
			Bins[i]->Fill(x,y);
		}else{
			OutOfBorder->Fill(x,y);
		}
	}
	
	Condition::Condition(ConditionTrackParamDependent func){
		condition=func;
	}
	Condition::Condition(ConditionIndependent func)
		:Condition([func](WTrack&,const vector<double>&){return func();}){}
	Condition::Condition(ConditionParamDependent func)
		:Condition([func](WTrack&,const vector<double>&P){return func(P);}){}
	Condition::Condition(ConditionTrackDependent func)
		:Condition([func](WTrack&T,const vector<double>&){return func(T);}){}
	Condition::~Condition(){}
	bool Condition::Process(WTrack&T, vector<double>&P)const{return condition(T,P);}
	
	Parameter::Parameter(ValueTrackParamDependent f){func=f;}
	Parameter::Parameter(ValueTrackDependent f):Parameter([f](WTrack&T,const vector<double>&){return f(T);}){}
	Parameter::Parameter(ValueParamDependent f):Parameter([f](WTrack&,const vector<double>&P){return f(P);}){}
	Parameter::Parameter(ValueIndependent f):Parameter([f](WTrack&,const vector<double>&){return f();}){}
	Parameter::~Parameter(){}
	bool Parameter::Process(WTrack&T, vector< double >&P)const{
		double x=func(T,P);
		P.push_back(x);
		return true;
	}
	
	AbstractChain::AbstractChain(){}
	AbstractChain::~AbstractChain(){}
	AbstractChain& AbstractChain::operator<<(shared_ptr< ITrackParamProcess > element){
		m_chain.push_back(element);
		return *this;
	}
	void AbstractChain::Cycle(function<void(bool)>func,WTrack&T,vector<double>&P) const{
		for(auto element:m_chain)
			func(element->Process(T,P));
	}
	bool ChainAnd::Process(WTrack&T, vector< double >&P) const{
		bool res=true;
		Cycle([&res](bool v){res&=v;},T,P);
		return res;
	}
	bool ChainOr::Process(WTrack&T, vector< double >&P) const{
		bool res=false;
		Cycle([&res](bool v){res|=v;},T,P);
		return res;
	}
	TrackProcess& TrackProcess::operator<<(shared_ptr< ITrackParamProcess > element){
		m_proc.push_back(element);
		return *this;
	}
	void TrackProcess::Process(WTrack& T) const{
		vector<double> P;
		for(auto proc:m_proc)
			proc->Process(T,P);
	}
}