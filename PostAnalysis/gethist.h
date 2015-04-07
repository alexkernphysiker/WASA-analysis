#ifndef DAIHADDG
#define DAIHADDG
#include <string>
#include <vector>
class hist{
public:
	struct point{
		double x,y,dx,dy;
	};
	hist(std::string filename,std::vector<std::string> &path,std::string histname);
	virtual ~hist();
	double Entries();
	point &operator[](int i);
	int count();
	typedef std::vector<point>::iterator iterator;
	typedef std::vector<point>::const_iterator const_iterator;
	iterator begin();
	const_iterator cbegin()const;
	iterator end();
	const_iterator cend() const;
private:
	std::vector<point> data;
	double norm;
};
#endif