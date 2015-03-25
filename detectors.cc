#include <vector>
#include <utility>
#include "detectors.h"
int LastForwardDepPlane(WTrack* track){
	static vector<pair<ForwardDetectorPlane,double>> planes;
	if(planes.size()==0){
		planes.push_back(make_pair(kFRH1,0.004));
		planes.push_back(make_pair(kFRH2,0.0025));
		planes.push_back(make_pair(kFRH3,0.0025));
		planes.push_back(make_pair(kFRH4,0.0035));
		planes.push_back(make_pair(kFRH5,0.004));
	}
	int res=0;int cnt=0;
	for(auto plane:planes){
		cnt++;
		if(track->Edep(plane.first)>plane.second){
			if(res<(cnt-1))
				return -1;//particle has tunneled through one of the planes :)
			cnt++;
		}
	}
	return res;
};

