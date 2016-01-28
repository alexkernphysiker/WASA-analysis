// this file is distributed under 
// MIT license
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
	Binner::Binner(double f, double t, unsigned int b){
		from=f;
		to=t;
		bins=b;
		CheckCorrectness();
	}
	Binner::Binner(const Binner& source):Binner(source.from,source.to,source.bins){}
	Binner& Binner::operator=(const Binner& source){
		from=source.from;
		to=source.to;
		bins=source.bins;
		CheckCorrectness();
		return *this;
	}
	Binner::~Binner(){}
	size_t Binner::count() const{return bins;}
	double Binner::left() const{return from;}
	double Binner::right() const{return to;}
	double Binner::bin_width() const{
		if(0==bins)throw Exception<Binner>("count==0");
		return (to-from)/double(bins);
	}
	double Binner::bin_center(size_t i) const{
		if(i>=bins)throw Exception<Binner>("bin range check error");
		return from+bin_width()*(double(i)+0.5);
	}
	void Binner::CheckCorrectness() const{
		if(to<=from)throw Exception<Binner>("wrong binning ranges");
		if(0==bins)throw Exception<Binner>("there cannot be zero bins");
	}

	bool Binner::FindBinIndex(size_t& output, double x) const{
		double delta=bin_width()/2.0;
		if(delta<=0)throw Exception<Binner>("delta<=0");
		for(size_t i=0,n=count();i<n;i++){
			double pos=bin_center(i);
			if((x>=(pos-delta))&&(x<(pos+delta))){
				output=i;
				return true;
			}
		}
		return false;
	}
	
	AHist1D::AHist1D(const string& dir, const string& name, const Binner& x){
		hist=new TH1F(name.c_str(),"",x.count(),x.left(),x.right());
		gHistoManager->Add(hist,dir.c_str());
	}
	AHist1D::~AHist1D(){}
	void AHist1D::Accept(double x)const{
		hist->Fill(x);
	}
	ASetOfHists1D::ASetOfHists1D(const string& dir, const string& name, const Binner& z, const Binner& x):Z(z){
		All=new TH1F((name+"-AllBins").c_str(),"",x.count(),x.left(),x.right());
		gHistoManager->Add(All,dir.c_str());
		OutOfBorder=new TH1F((name+"-OutOfBins").c_str(),"",x.count(),x.left(),x.right());
		gHistoManager->Add(OutOfBorder,dir.c_str());
		for(unsigned int i=0;i<Z.count();i++){
			TH1F*hist=new TH1F((name+"-Bin-"+to_string(i)).c_str(),"",x.count(),x.left(),x.right());
			gHistoManager->Add(hist,dir.c_str());
			Bins.push_back(hist);
		}
	}
	ASetOfHists1D::~ASetOfHists1D(){}
	void ASetOfHists1D::Accept(double z, double x) const{
		size_t i=0;
		if(Z.FindBinIndex(i,z)){
			All->Fill(x);
			Bins[i]->Fill(x);
		}else{
			OutOfBorder->Fill(x);
		}
	}
	AHist2D::AHist2D(const string&dir,const string& name, const Binner& x, const Binner& y){
		hist=new TH2F(name.c_str(),"",x.count(),x.left(),x.right(),y.count(),y.left(),y.right());
		gHistoManager->Add(hist,dir.c_str());
	}
	AHist2D::~AHist2D(){}
	void AHist2D::Accept(double x, double y) const{
		hist->Fill(x,y);
	}
	ASetOfHists2D::ASetOfHists2D(const string&dir,const string& name, const Binner& z, const Binner& x, const Binner& y):Z(z){
		All=new TH2F((name+"-AllBins").c_str(),"",x.count(),x.left(),x.right(),y.count(),y.left(),y.right());
		gHistoManager->Add(All,dir.c_str());
		OutOfBorder=new TH2F((name+"-OutOfBins").c_str(),"",x.count(),x.left(),x.right(),y.count(),y.left(),y.right());
		gHistoManager->Add(OutOfBorder,dir.c_str());
		for(unsigned int i=0;i<Z.count();i++){
			TH2F*hist=new TH2F((name+"-Bin-"+to_string(i)).c_str(),"",x.count(),x.left(),x.right(),y.count(),y.left(),y.right());
			gHistoManager->Add(hist,dir.c_str());
			Bins.push_back(hist);
		}
	}
	ASetOfHists2D::~ASetOfHists2D(){}
	void ASetOfHists2D::Accept(double z, double x, double y) const{
		size_t i=0;
		if(Z.FindBinIndex(i,z)){
			All->Fill(x,y);
			Bins[i]->Fill(x,y);
		}else{
			OutOfBorder->Fill(x,y);
		}
	}

	class ParamlessConverter:public ICalculation{
	private:
		shared_ptr<ParamlessProcess> func;
	public:
		ParamlessConverter(shared_ptr<ParamlessProcess> f){func=f;}
		virtual ~ParamlessConverter(){}
		virtual bool Process(WTrack&,Results&)const override{
			return func->Process();
		}
	};
	class AnalysisConverter:public ICalculation{
	private:
		shared_ptr<IResAnalysis> func;
	public:
		AnalysisConverter(shared_ptr<IResAnalysis> f){func=f;}
		virtual ~AnalysisConverter(){}
		virtual bool Process(WTrack&,Results&P)const override{
			return func->Process(P);
		}
	};
	class TrackConverter:public ICalculation{
	private:
		shared_ptr<ITrackProcess> func;
	public:
		TrackConverter(shared_ptr<ITrackProcess> f){func=f;}
		virtual ~TrackConverter(){}
		virtual bool Process(WTrack&T,Results&)const override{
			return func->Process(T);
		}
	};
	void CalculationSet::AddC(shared_ptr<ParamlessProcess>s){
		AbstractSet<WTrack&,Results&>::Add(make_shared<ParamlessConverter>(s));
	}
	void CalculationSet::AddC(shared_ptr<IResAnalysis> s){
		AbstractSet<WTrack&,Results&>::Add(make_shared<AnalysisConverter>(s));
	}
	void CalculationSet::AddC(shared_ptr<ITrackProcess>s){
		AbstractSet<WTrack&,Results&>::Add(make_shared<TrackConverter>(s));
	}
	bool TrackCalculations::Process(WTrack& T) const{
		Results P;
		for(auto p:allset())p->Process(T,P);
		return true;
	}
	bool TrackCalcChain::Process(WTrack& T) const{
		Results P;
		for(auto p:allset())p->Process(T,P);
		return true;
	}


}