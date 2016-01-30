// this file is distributed under 
// MIT license
#ifndef HGDOYENRKILTMEYD
#define HGDOYENRKILTMEYD
#include <string>
#include <functional>
#include <TH2F.h>
#include "log.h"
#include "trackprocessing.h"
enum ForwardDetectorPlane{
	kForwardError=0,
	kFWC1 = 10, kFWC2 = 11, kFTH1 = 1, kFTH2 =2, kFTH3 = 3, 
	kFRH1 = 4, kFRH2 = 5, kFRH3 = 6, kFRH4 = 7,kFRH5 = 8, kFVH =9
};
namespace TrackAnalyse{
	using namespace std;
	class Forward{
	protected:
		Forward();
	public:
		static const Forward&Get();
		virtual ~Forward();
		class plane_data{
			friend class Forward;
		public:
			plane_data(ForwardDetectorPlane p,string&&n,double thr,double upper);
			~plane_data();
			double Edep(WTrack&)const;
		protected:
			ForwardDetectorPlane plane()const;
			string&name()const;
			Axis&axis()const;
			double threshold()const;
		private:
			ForwardDetectorPlane m_plane;
			string m_name;
			double m_thr;
			Axis m_axis;
		};
		typedef vector<plane_data>::const_iterator const_iterator;
		const_iterator begin()const;
		const_iterator cbegin()const;
		const_iterator end()const;
		const_iterator cend()const;
		size_t count()const;
		plane_data&operator[](ForwardDetectorPlane)const;
		ForwardDetectorPlane StoppingLayer(WTrack&)const;
		shared_ptr<Chain> CreateMarker(string&&dir,string&&name)const;
	private:
		vector<plane_data> PlaneData;
	};
}
#endif
