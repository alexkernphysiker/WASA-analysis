// this file is distributed under 
// MIT license
#include "math_h/error.h"
#include "trackprocessing.h"
#include "../config.h"
using namespace std;
using namespace MathTemplates;
TrackConditionSet::TrackConditionSet(string&&name, ValueIndependent distr,int bins,double from,double to){
	m_name=name;m_distr=distr;m_bins=bins;m_from=from;m_to=to;
	reference=new TH1F("Reference","",m_bins,m_from,m_to);
	gHistoManager->Add(reference,m_name.c_str());
	beforecut=new TH1F("BeforeAllCuts","",m_bins,m_from,m_to);
	gHistoManager->Add(beforecut,m_name.c_str());
}
TrackConditionSet::TrackConditionSet(string&&name, const TrackConditionSet& master):
	TrackConditionSet(static_cast<string&&>(name),master.m_distr,master.m_bins,master.m_from,master.m_to){}
TrackConditionSet::~TrackConditionSet(){}
TrackConditionSet::TrackCalc::TrackCalc(const string&n):name(n){}
TrackConditionSet::TrackCalc::~TrackCalc(){}
string&& TrackConditionSet::TrackCalc::Name(){return static_cast<string&&>(name);}
TrackConditionSet::ParamCalc::ParamCalc(const string&n, ValueTrackDependent d):TrackCalc(n){m_delegate=d;}
TrackConditionSet::ParamCalc::~ParamCalc(){}
double TrackConditionSet::ParamCalc::Get(WTrack&track){
	return m_delegate(track);
}
TrackConditionSet::TrackCondition::TrackCondition(const string&n,InternalCondition delegate, TrackConditionSet* master):
	TrackCalc(n){
	condition=delegate;
	output=new TH1F(n.c_str(),"",master->m_bins,master->m_from,master->m_to);
	gHistoManager->Add(output,master->m_name.c_str());
}
TrackConditionSet::TrackCondition::~TrackCondition(){}
bool TrackConditionSet::TrackCondition::Check(WTrack&track,vector<double>&P,double magnitude){
	if(!condition(track,P))
		return false;
	output->Fill(magnitude);
	return true;
}
TrackConditionSet& TrackConditionSet::AddParameter(std::string&&n,ValueTrackDependent parameter){
	calc_procs.push_back(make_shared<ParamCalc>(n,parameter));
	return *this;
}
TrackConditionSet& TrackConditionSet::AddCondition(string&&n,ConditionTrackParamDependent condition){
	calc_procs.push_back(make_shared<TrackCondition>(n,[condition](WTrack&T,vector<double>&P){return condition(T,P);},this));
	return *this;
}
TrackConditionSet& TrackConditionSet::AddCondition(string&& name, ConditionParamDependent condition){
	return AddCondition(static_cast<string&&>(name),[condition](WTrack&,const vector<double>&P){return condition(P);});
}
TrackConditionSet& TrackConditionSet::AddCondition(string&& name, ConditionTrackDependent condition){
	return AddCondition(static_cast<string&&>(name),[condition](WTrack&T,const vector<double>&){return condition(T);});
}
TrackConditionSet& TrackConditionSet::AddCondition(string&& name, ConditionIndependent condition){
	return AddCondition(static_cast<string&&>(name),[condition](WTrack&,const vector<double>&){return condition();});
}
TrackConditionSet& TrackConditionSet::AddConditions(string&& name,vector< TrackConditionSet >& set){
	calc_procs.push_back(make_shared<TrackCondition>(name,[&set](WTrack&T,vector<double>&P){
		bool res=false;
		for(TrackConditionSet&cond:set)res|=cond.Check(T,P);
		return res;
	},this));
	return *this;
}
void TrackConditionSet::ReferenceEvent(){
	double M=m_distr();
	reference->Fill(M);
}
bool TrackConditionSet::Check(WTrack&track,vector<double>&params){
	double M=m_distr();
	beforecut->Fill(M);
	for(auto item:calc_procs){
		if(dynamic_pointer_cast<ParamCalc>(item))
			 params.push_back(dynamic_pointer_cast<ParamCalc>(item)->Get(track));
		if(dynamic_pointer_cast<TrackCondition>(item))
			if(!dynamic_pointer_cast<TrackCondition>(item)->Check(track,params,M))
				return false;
	}
	return true;
}

Analyser2D::Analyser2D(std::string&&name,const TrackConditionSet&Cuts){
	m_name=name;
	ForAllA=nullptr;
	reference=Cuts.reference;
	m_distr=Cuts.m_distr;
}
void Analyser2D::Setup(ValueParamDependent B, int binsB, double fromB, double toB){
	m_B=B;
	ForAllA=new TH1F("All","",binsB,fromB,toB);
	gHistoManager->Add(ForAllA,m_name.c_str());
	{
		TH1F* OutOff=new TH1F("Out","",binsB,fromB,toB);
		A_bin.push_back(OutOff);
		gHistoManager->Add(OutOff,m_name.c_str());
	}
	for(int i=1,N=reference->GetNbinsX();i<=N;i++){
		int index=int(reference->GetBinCenter(i)*binning_coefficient);
		TH1F* hist=new TH1F(to_string(index).c_str(),"",binsB,fromB,toB);
		A_bin.push_back(hist);
		gHistoManager->Add(hist,m_name.c_str());
	}
}
Analyser2D::~Analyser2D(){}
void Analyser2D::AcceptEvent(const vector<double>&Parameters){
	double B=m_B(Parameters);
	double A=m_distr();
	ForAllA->Fill(B);
	{int index=0;
		for(int i=1,N=reference->GetNbinsX();(i<=N)&&(index==0);i++){
			double min=reference->GetBinLowEdge(i);
			if((A>=min)&&(A<(min+reference->GetBinWidth(i))))
				index=i;
		}
		A_bin[index]->Fill(B);
	}
}

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

Debug2DSpectraSet::Debug2DSpectraSet(string&&name):m_name(name){}
Debug2DSpectraSet::~Debug2DSpectraSet(){}
Debug2DSpectraSet::Process Debug2DSpectraSet::Create(ValueTrackParamDependent x, ValueTrackParamDependent y){
	return [this,x,y](WTrack&track,const vector<double>&P){return make_pair(x(track,P),y(track,P));};
}
void Debug2DSpectraSet::CatchState(WTrack& track,const vector<double>&params){
	for(Item&item:jobs){
		point P=item.second(track,params);
		item.first->Fill(P.first,P.second);
	}
}
void Debug2DSpectraSet::Add(string&&name,const Axis& X, const Axis& Y){
	TH2F* hist=new TH2F(name.c_str(),"",X.count(),X.left(),X.right(),Y.count(),Y.left(),Y.right());
	gHistoManager->Add(hist,m_name.c_str());
	jobs.push_back(make_pair(hist,Create(X.valuegetter(),Y.valuegetter())));
}
void Debug2DSpectraSet::Add(string&& name, Axis&& X, Axis&& Y){
	Add(static_cast<string&&>(name),X,Y);
}

Debug2DSpectraBined::Debug2DSpectraBined(string&& name, const Axis& x, const Axis& y, const Axis& z):X(x),Y(y),Z(z){
	for(unsigned int i=0; i<Z.count();i++){
		TH2F* hist=new TH2F(to_string(Z.bin_center(i)).c_str(),"",X.count(),X.left(),X.right(),Y.count(),Y.left(),Y.right());
		gHistoManager->Add(hist,name.c_str());
		sp2.push_back(hist);
	}
}
Debug2DSpectraBined::Debug2DSpectraBined(string&& name, Axis&& x, Axis&& y, Axis&& z):Debug2DSpectraBined(static_cast<string&&>(name),x,y,z){}
Debug2DSpectraBined::~Debug2DSpectraBined(){}
void Debug2DSpectraBined::CatchState(WTrack& T, const vector< double >& P){
	unsigned int i=0;
	if(Z.FindBinIndex(i,T,P))
		sp2[i]->Fill(X.getvalue(T,P),Y.getvalue(T,P));
}
