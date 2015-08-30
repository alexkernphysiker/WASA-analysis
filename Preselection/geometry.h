// this file is distributed under 
// GPL v 3.0 license
#ifndef gYQHDBYV
#define gYQHDBYV
#include <exception>
#include "math_h/interpolate.h"
namespace PlaneGeometry{
#define Point pair<numX,numX>
#define Segment pair<pair<numX,numX>,pair<numX,numX>>
	using namespace std;
	template<class numX>
	inline numX CrossProduct(const Point&A,const Point&B){
		return B.first*A.second-A.first*B.second;
	}
	template<class numX>
	inline Point Vector(const Segment&S){
		return make_pair(S.first.first-S.second.first,S.first.second-S.second.second);
	}
	
	template<class numX>
	bool SegmentsIntersect(const Segment&A,const Segment&B){
		auto Va=Vector(A);
		if(CrossProduct<numX>(Va,B.first)*CrossProduct<numX>(Va,B.second)<=0){
			auto Vb=Vector(B);
			if(CrossProduct<numX>(Vb,A.first)*CrossProduct<numX>(Vb,A.second)<=0)
				return true;
			else
				return false;
		}else 
			return false;
	}
	
	template<class numX>class Polygon{
	private:
		vector<Point> data;
	public:
		Polygon(){}
		virtual Polygon &operator<<(Point&&p){
			size_t n=data.size();
			if(n>1){
				Segment last1=make_pair(p,data[0]);
				Segment last2=make_pair(data[n-1],p);
				for(size_t i=1;i<n;i++){
					auto cur=make_pair(data[i-1],data[i]);
					if(SegmentsIntersect<numX>(last1,cur))
						throw exception();
					if(SegmentsIntersect<numX>(last2,cur))
						throw exception();
				}
			}
			data.push_back(p);
			return *this;
		}
		Polygon(vector<Point>&&points){
			for(auto&p:points)
				operator<<(static_cast<Point&&>(p));
		}
		Polygon(const Polygon&source){
			for(auto&p:source.data)
				operator<<(static_cast<Point&&>(p));
		}
		virtual ~Polygon(){}
		size_t size()const{return data.size();}
		Point&operator[](int i)const{return data[i];}
		typedef typename vector<Point>::iterator iterator;
		typedef typename vector<Point>::const_iterator const_iterator;
		iterator begin(){return data.begin();}
		const_iterator begin()const{return data.begin();}
		const_iterator cbegin()const{return data.cbegin();}
		iterator end(){return data.end();}
		const_iterator end() const{return data.end();}
		const_iterator cend() const{return data.cend();}
		
		bool operator()(Point&&X){
			//ToDo: implement
			return false;
		}
	};
#undef Point
#undef Segment
};
#endif 