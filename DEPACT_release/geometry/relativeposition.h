/*
 * relativeposition.h
 *
 *  Created on: 2016年4月1日
 *      Author: hyliu
 */

#ifndef RELATIVEPOSITION_H_
#define RELATIVEPOSITION_H_
#include "geometry/xyz.h"
#include "geometry/quaternioncrd.h"
#include <iostream>
namespace NSPgeometry{
struct RelativePosition {
	XYZ location;
	QuaternionCrd orientation;
	void read(std::istream & is) {
		is >>location.x_ >> location.y_ >> location.z_
		   >> orientation[0] >> orientation[1] >>orientation[2] >>orientation[3];
	}
	void write(std::ostream &os) {
		os <<" " <<location.x_ <<" " <<location.y_<<" " <<location.z_ <<" "
				<< orientation[0]<<" " <<orientation[1]<<" " << orientation[2]
				<<" " <<orientation[3]<< std::endl;
	}
	std::vector<double> tovector() {
		std::vector<double> vec=location.tovector();
		for(int d=0; d<4;d++) {
			vec.push_back(orientation[d]);
		}
	}
};

template <typename RNG>
RelativePosition fixdis_randomorientation (double r, RNG &rng) {
	RelativePosition rp;
	bool get=false;
	while(!get) {
		double norm=0.0;
		std::vector<double> xyz;
		for(int d=0; d<3; d++) {
			double x=rng(); //random uniform distribution between 0 and 1.
			x= x*2.0-1.0;
			xyz.push_back(x);
			norm +=x*x;
		}
		if(norm >1.0) continue;
		get=true;
		XYZ p(xyz);
		rp.location=(r/sqrt(p.squarednorm()))*p;
	}
	get=false;
	while (!get) {
		double norm=0.0;
		std::vector<double> q;
		for(int d=0; d<4; ++d) {
			double x=rng();
			x=x*2.0-1.0;
			q.push_back(x);
			norm +=x*x;
		}
		if(norm >1.0) continue;
		get=true;
		norm=1.0/sqrt(norm);
		for(int d=0; d<4; ++d) q[d] *=norm;
		rp.orientation=QuaternionCrd(q);
	}
	return rp;
}

template <typename RNG>
RelativePosition samplerelativeposition (double rmax, RNG &rng) {
	RelativePosition rp;
	bool get=false;
	while(!get) {
		double norm=0.0;
		std::vector<double> xyz;
		for(int d=0; d<3; d++) {
			double x=rng(); //random uniform distribution between 0 and 1.
			x= x*2.0-1.0;
			xyz.push_back(x);
			norm +=x*x;
		}
		if(norm >1.0) continue;
		get=true;
		rp.location=rmax*XYZ(xyz);
	}
	get=false;
	while (!get) {
		double norm=0.0;
		std::vector<double> q;
		for(int d=0; d<4; ++d) {
			double x=rng();
			x=x*2.0-1.0;
			q.push_back(x);
			norm +=x*x;
		}
		if(norm >1.0) continue;
			get=true;
			norm=1.0/sqrt(norm);
			for(int d=0; d<4; ++d) q[d] *=norm;
			rp.orientation=QuaternionCrd(q);
	}
	return rp;
}
}

#endif /* RELATIVEPOSITION_H_ */
