// this file is distributed under 
// MIT license
#include "trackprocessing.h"
#include "../config.h"
using namespace std;
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
	value=v;
	from=f;
	to=t;
	bins=b;
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
	TH2F* hist=new TH2F(name.c_str(),"",X.bins,X.from,X.to,Y.bins,Y.from,Y.to);
	gHistoManager->Add(hist,m_name.c_str());
	jobs.push_back(make_pair(hist,Create(X.value,Y.value)));
}
