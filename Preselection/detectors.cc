// this file is distributed under 
// MIT license
#include <vector>
#include <utility>
#include "math_h/error.h"
#include "detectors.h"
namespace TrackAnalyse {
	using namespace std;
	ForwardDetectors::plane_data::plane_data(ForwardDetectorPlane p, string&& n, double thr,double upper)
		:m_plane(p),m_name(n),m_thr(thr),m_axis([p](WTrack&T){return T.Edep(p);},0,upper,100){}
	ForwardDetectors::plane_data::~plane_data(){}
	ForwardDetectorPlane ForwardDetectors::plane_data::plane() const{return m_plane;}
	Axis& ForwardDetectors::plane_data::axis() const{return const_cast<Axis&>(m_axis);}
	string& ForwardDetectors::plane_data::name() const{return const_cast<string&>(m_name);}
	double ForwardDetectors::plane_data::threshold() const{return m_thr;}
	double ForwardDetectors::plane_data::Edep(WTrack&T) const{
		vector<double> P;
		return m_axis.getvalue(T,P);
	}
	
	ForwardDetectors::ForwardDetectors(){
		PlaneData.push_back(plane_data(kFWC1,"FWC1",0.002 ,0.03));
		PlaneData.push_back(plane_data(kFWC2,"FWC2",0.002 ,0.03));
		//PlaneData.push_back(plane_data(kFPC,"FPC",0.002 ,0.03));
		PlaneData.push_back(plane_data(kFTH1,"FTH1",0.0015,0.05));
		PlaneData.push_back(plane_data(kFRH1,"FRH1",0.001, 0.05));
		PlaneData.push_back(plane_data(kFRH2,"FRH2",0.001, 0.05));
	}
	ForwardDetectors::~ForwardDetectors(){}
	const ForwardDetectors& ForwardDetectors::Instance(){
		static ForwardDetectors result;
		return const_cast<ForwardDetectors&>(result);
	}
	ForwardDetectors::const_iterator ForwardDetectors::begin() const{
		return PlaneData.begin();
	}
	ForwardDetectors::const_iterator ForwardDetectors::cbegin() const{
		return PlaneData.cbegin();
	}
	ForwardDetectors::const_iterator ForwardDetectors::end() const{
		return PlaneData.end();
	}
	ForwardDetectors::const_iterator ForwardDetectors::cend() const{
		return PlaneData.cend();
	}
	size_t ForwardDetectors::count() const{return PlaneData.size();}
	ForwardDetectors::plane_data& ForwardDetectors::operator[](ForwardDetectorPlane index) const{
		for(const plane_data&plane:PlaneData)
			if(plane.plane()==index)
				return const_cast<plane_data&>(plane);
		throw MathTemplates::Exception<ForwardDetectors>("Forward layer not found");
	}
	ForwardDetectorPlane ForwardDetectors::StoppingLayer(WTrack&T)const{
		int res=-1;
		for(size_t i=0;i<PlaneData.size();i++)
			if(PlaneData[i].Edep(T)>PlaneData[i].threshold()){
				if(res==(int(i)-1))res++;
				else res=-1;
			}
		if(res>=0)
			return PlaneData[res].plane();
		return kForwardError;
	}
	shared_ptr<Chain> ForwardDetectors::CreateMarker(string&& dir, string&& name)const{
		auto res=make_shared<Chain>();
		for(size_t i=1;i<PlaneData.size();i++)
			res << make_shared<Hist2D>(
				static_cast<string&&>(dir),name+"-"+PlaneData[i-1].name()+"-vs-"+PlaneData[i].name(),
				PlaneData[i].axis(),PlaneData[i-1].axis()
			);
		return res;
	}
}
