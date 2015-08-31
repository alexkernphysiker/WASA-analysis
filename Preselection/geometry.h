// this file is distributed under 
// GPL v 3.0 license
#ifndef gYQHDBYV
#define gYQHDBYV
#include <exception>
#include "math_h/interpolate.h"
namespace PlaneGeometry{
#define Point pair<numX,numX>
#define Segment pair<pair<numX,numX>,pair<numX,numX>>
	template<class numX>class exception:public std::exception{
	private:
		std::string m_msg;
	public:
		exception<numX>(std::string&&msg){m_msg=msg;}
		virtual ~exception<numX>() throw(){}
		virtual const char* what() const throw(){return m_msg.c_str();}
	};
	using namespace std;
	template<class numX>inline bool operator==(const Point&A,const Point&B){
		return (A.first==B.first)&&(A.second==B.second);
	}
	template<class numX>inline bool operator!=(const Point&A,const Point&B){
		return (A.first!=B.first)||(A.second!=B.second);
	}
	template<class numX>inline numX CrossProduct(const Point&A,const Point&B){
		return B.first*A.second-A.first*B.second;
	}
	template<class numX>inline numX CrossProduct(Point&&A,Point&&B){
		return CrossProduct<numX>(const_cast<Point&>(A),const_cast<Point&>(B));
		
	}
	template<class numX>inline Point operator+(const Point&A,const Point&B){
		return make_pair<numX,numX>(A.first+B.first,A.second+B.second);
	}
	template<class numX>inline Point operator+(Point&&A,Point&&B){
		return make_pair<numX,numX>(A.first+B.first,A.second+B.second);
	}
	template<class numX>inline Point operator-(const Point&A,const Point&B){
		return make_pair<numX,numX>(A.first-B.first,A.second-B.second);
	}
	template<class numX>inline Point operator-(Point&&A,Point&&B){
		return make_pair<numX,numX>(A.first-B.first,A.second-B.second);
	}
	template<class numX>inline Point operator*(const Point&A,numX B){
		return make_pair<numX,numX>(A.first*B,A.second*B);
	}
	template<class numX>inline Point operator*(Point&&A,numX B){
		return make_pair<numX,numX>(A.first*B,A.second*B);
	}
	template<class numX>inline numX operator*(const Point&A,const Point&B){
		return A.first*B.first+A.second*B.second;
	}
	template<class numX>inline numX operator*(Point&&A,Point&&B){
		return A.first*B.first+A.second*B.second;
	}
	
	template<class numX>inline bool PointOnLine(const Point&p,const Segment&S){
		return CrossProduct<numX>(p-S.first,S.second-S.first)==0;
	}
	template<class numX>
	bool SegmentsIntersect(const Segment&A,const Segment&B){
		if(CrossProduct(A.second-A.first,B.first-A.first)*CrossProduct<numX>(A.second-A.first,B.second-A.first)<=0)
			if(CrossProduct(B.second-B.first,A.first-B.first)*CrossProduct<numX>(B.second-B.first,A.second-B.first)<=0)
				return true;
		return false;
	}
	template<class numX>
	bool SegmentIntersectsLine(const Segment&S,const Segment&L){
		return (CrossProduct(L.second-L.first,S.first-L.first)*CrossProduct<numX>(L.second-L.first,S.second-L.first)<=0);
	}
	template<class numX>numX LineYbyX(const Segment&A,double X){
		if(A.first.first==A.second.first) throw exception<numX>("LineYbyX error");
		return (X-A.first.first)*(A.second.second-A.first.second)/(A.second.first-A.first.first);
	}
	template<class numX>numX LineXbyY(const Segment&A,double Y){
		if(A.first.second==A.second.second) throw exception<numX>("LineXbyY error");
		return (Y-A.first.second)*(A.second.first-A.first.first)/(A.second.second-A.first.second);
	}
	
	template<class numX>class Polygon{
	private:
		vector<Point> data;
	public:
		Polygon(){}
		Polygon &operator<<(const Point&p){
			size_t n=data.size();
			if(n>1){
				Segment last=make_pair(data[n-1],p);
				for(size_t i=1;i<n;i++){
					Segment cur=make_pair(data[i-1],data[i]);
					if((i>1)||(p!=data[0]))if(SegmentsIntersect<numX>(last,cur))
						throw exception<numX>("Polygon: intersection detected!");
				}
			}
			data.push_back(p);
			return *this;
		}
		Polygon &operator<<(Point&&p){return operator<<(const_cast<Point&>(p));}
		Polygon(const Polygon&source){
			for(auto&p:source.data)operator<<(p);
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
			vector<numX> online;
			size_t n=size();
			if(data[0]!=data[n-1])
				throw exception<numX>("Using of incomplete polygon");
			for(size_t i=1;i<n;i++){
				Segment seg=make_pair(data[i-1],data[i]);
				if(seg.second.first!=seg.first.first)
					if((seg.first.first-X.first)*(seg.second.first-X.first)<=0)
						InsertSorted(LineYbyX(seg,X.first),online,std_size(online),std_insert(online,numX));
			}
			if(online.size()%2==1)throw exception<numX>("Polynom: Odd number of intersections with vertical line");
			size_t index=0;
			while((index<online.size())&&(X.second<=online[index]))index++;
			return (index%2)==1;
		}
	};
#undef Point
#undef Segment
};
#endif 