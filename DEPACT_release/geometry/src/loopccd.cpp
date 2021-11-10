/*
 * loopccd.cpp
 *
 *  Created on: 2016年11月23日
 *      Author: hyliu
 */

#include "geometry/loopccd.h"
#include "geometry/localframe.h"
#include "geometry/rotation.h"
#include <cassert>
using namespace NSPgeometry;

double NSPgeometry::ccdangle(const Line &axis, const std::vector<XYZ> &fxpoints, const std::vector<XYZ> &mvpoints){
	const double EPSINON=1.e-20;
	const double PI=3.14159265358979323846;
	assert(fxpoints.size() == mvpoints.size());
	double a=0.0,b=0.0;
	for (int i=0; i< fxpoints.size(); ++i) {
		LocalFrame lf=make_localframe(axis,mvpoints[i]);
		XYZ p1=lf.global2localcrd(fxpoints[i]);
		double rxy=sqrt(p1.x_*p1.x_+p1.y_*p1.y_);
		a += p1.y_*rxy;
		b +=p1.x_*rxy;
	}
	return atan2(a,b);
}


double NSPgeometry::rotationtormsd(const Line &axis, double rotation,const std::vector<XYZ> &fxpoints, std::vector<XYZ> mvpoints){
	Rotation R;
	const double PI=3.14159265358979323846;
	R.init(QuaternionCrd(axis.direction, rotation*180.0/PI),axis.origin);
	assert(fxpoints.size() == mvpoints.size());
	double newrmsd=0.0;
	for(auto &r:mvpoints) R.apply(&r);
	newrmsd= NSPgeometry::rmsd(fxpoints,mvpoints);
	return newrmsd;
}

