// this file is distributed under 
// GPL v 3.0 license
#ifndef gYQHDBYV
#define gYQHDBYV
#include "math_h/interpolate.h"
template<class numX,class numY=numX>
class Polygon{
public:
	typedef std::pair<numX,numY> Point;
private:
	std::vector<Point> data;
public:
	Polygon(){}
	virtual Polygon &operator<<(Point&&p){
		data.push_back(p);
		return *this;
	}
	Polygon(std::vector<Point>&&points){
		for(Point&p:points)
			operator<<(static_cast<Point&&>(p));
	}
	Polygon(const Polygon&source){
		for(Point&p:source.data)
			operator<<(static_cast<Point&&>(p));
	}
	virtual ~Polygon(){}
	int size()const{return data.size();}
	Point&operator[](int i)const{return data[i];}
	typedef typename std::vector<Point>::iterator iterator;
	typedef typename std::vector<Point>::const_iterator const_iterator;
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
	bool operator()(numX x,numY y){return operator()(std::make_pair(x,y));}
};
#endif 