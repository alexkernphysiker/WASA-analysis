// this file is distributed under 
// GPL v 3.0 license
#include "trackprocessing.h"
#define static_right(A) (static_cast<decltype(A)&&>(A))
using namespace std;
TrackConditionSet::TrackConditionSet(string name, Independent distr,int bins,double from,double to){
	m_name=name;m_distr=distr;
	m_from=from;m_to=to;m_bins=bins;
	reference=new TH1F("Reference","",m_bins,m_from,m_to);
	gHistoManager->Add(reference,m_name.c_str());
	beforecut=new TH1F("BeforeAllCuts","",m_bins,m_from,m_to);
	gHistoManager->Add(beforecut,m_name.c_str());
}
TrackConditionSet::~TrackConditionSet(){}
TrackConditionSet& TrackConditionSet::AddParameter(TrackDependent parameter){
	parameter_set.push_back(parameter);
	return *this;
}
TrackConditionSet::CondData::CondData(string n, TrackConditionSet::Condition delegate, TrackConditionSet* master){
	condition=delegate;
	output=new TH1F(n.c_str(),"",master->m_bins,master->m_from,master->m_to);
	gHistoManager->Add(output,master->m_name.c_str());
}
TrackConditionSet& TrackConditionSet::AddCondition(string name, TrackConditionSet::Condition condition){
	condition_set.push_back(CondData(name,condition,this));
	return *this;
}
void TrackConditionSet::ReferenceEvent(){
	double M=m_distr();
	reference->Fill(M);
}
bool TrackConditionSet::Check(WTrack&& track,vector<double>&params){
	params.clear();
	for(TrackDependent par:parameter_set)
		params.push_back(par(static_right(track)));
	double M=m_distr();
	beforecut->Fill(M);
	for(CondData&item:condition_set){
		if(!item.condition(static_right(track),params))
			return false;
		item.output->Fill(M);
	}
	return true;
}

Analyser2D::Analyser2D(std::string name,TrackConditionSet&&Cuts,ParamDependent B,int binsB,double fromB,double toB){
	ForAllA=new TH1F("All","",binsB,fromB,toB);
	gHistoManager->Add(ForAllA,name.c_str());
	m_B=B;
	master=&Cuts;
	{
		TH1F* OutOff=new TH1F("Out","",binsB,fromB,toB);
		A_bin.push_back(OutOff);
		gHistoManager->Add(OutOff,name.c_str());
	}
	for(int i=1,N=Cuts.reference->GetNbinsX();i<=N;i++){
		int index=int(Cuts.reference->GetBinCenter(i)*1000);
		TH1F* hist=new TH1F(to_string(index).c_str(),"",binsB,fromB,toB);
		A_bin.push_back(hist);
		gHistoManager->Add(hist,name.c_str());
	}
}
Analyser2D::~Analyser2D(){}
void Analyser2D::AcceptEvent(vector<double>&&Parameters){
	double B=m_B(static_right(Parameters));
	double A=master->m_distr();
	ForAllA->Fill(B);
	{int index=0;
		for(int i=1,N=master->reference->GetNbinsX();(i<=N)&&(index==0);i++){
			double min=master->reference->GetBinLowEdge(i);
			if((A>=min)&&(A<(min+master->reference->GetBinWidth(i))))
				index=i;
		}
		A_bin[index]->Fill(B);
	}
}
