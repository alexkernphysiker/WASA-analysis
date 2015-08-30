// this file is distributed under 
// GPL v 3.0 license
#include "trackprocessing.h"
using namespace std;
TrackConditionSet::TrackConditionSet(string&&name, Independent distr,int bins,double from,double to){
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
TrackConditionSet::ParamCalc::ParamCalc(const string&n, TrackDependent d):TrackCalc(n){m_delegate=d;}
TrackConditionSet::ParamCalc::~ParamCalc(){}
double TrackConditionSet::ParamCalc::Get(WTrack&track){
	return m_delegate(track);
}
TrackConditionSet::TrackCondition::TrackCondition(const string&n,TrackConditionSet::Condition delegate, TrackConditionSet* master):
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
TrackConditionSet& TrackConditionSet::AddParameter(std::string&&n,TrackDependent parameter){
	calc_procs.push_back(make_shared<ParamCalc>(n,parameter));
	return *this;
}
TrackConditionSet& TrackConditionSet::AddCondition(string&&n, TrackConditionSet::Condition condition){
	calc_procs.push_back(make_shared<TrackCondition>(n,condition,this));
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
void Analyser2D::Setup(ParamDependent B, int binsB, double fromB, double toB){
	m_B=B;
	ForAllA=new TH1F("All","",binsB,fromB,toB);
	gHistoManager->Add(ForAllA,m_name.c_str());
	{
		TH1F* OutOff=new TH1F("Out","",binsB,fromB,toB);
		A_bin.push_back(OutOff);
		gHistoManager->Add(OutOff,m_name.c_str());
	}
	for(int i=1,N=reference->GetNbinsX();i<=N;i++){
		int index=int(reference->GetBinCenter(i)*1000);
		TH1F* hist=new TH1F(to_string(index).c_str(),"",binsB,fromB,toB);
		A_bin.push_back(hist);
		gHistoManager->Add(hist,m_name.c_str());
	}
}
Analyser2D::~Analyser2D(){}
void Analyser2D::AcceptEvent(vector<double>&Parameters){
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
