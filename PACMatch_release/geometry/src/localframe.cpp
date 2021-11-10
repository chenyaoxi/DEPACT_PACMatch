/*
 * localframe.cpp
 *
 *  Created on: 2016年3月31日
 *      Author: hyliu
 */

#include "geometry/localframe.h"
#include "geometry/quaternioncrd.h"
#include <cmath>

using namespace NSPgeometry;
A3LocalFrame::A3LocalFrame (XYZ a1, XYZ a2,XYZ a3){
	LocalFrame lf=make_localframe(a1,a2,a3);
	origin_=lf.origin_;
	axis_=lf.axis_;
	dlfda.assign(3,DLocalFrameDx());
	std::vector<XYZ> a123;
	a123.push_back(a1);
	a123.push_back(a2);
	a123.push_back(a3);
	for(int a=0;a<3;++a){
		dlfda[a].daxisdx.assign(3,std::vector<XYZ>());
		for(int m=0;m<3;++m){
			a123[a][m] +=0.00005;
			LocalFrame lfp=make_localframe(a123[0],a123[1],a123[2]);
			a123[a][m] -=0.0001;
			LocalFrame lfm=make_localframe(a123[0],a123[1],a123[2]);
			dlfda[a].doridx.push_back((lfp.origin_-lfm.origin_)/0.0001);
			for(int k=0;k<3;++k){
				dlfda[a].daxisdx[k].push_back((lfp.axis_[k]-lfm.axis_[k])/0.0001);
			}
			a123[a][m] +=0.00005;
		}
	}
}
std::vector<XYZ> A3LocalFrame::distributedvdx(const XYZ  & x, const XYZ & dvdx) const {
	XYZ localx=global2localcrd(x);
//	x=origin_+localx[0]*axis_[0]+localx[1]*axis_[1]+localx[2]*axis_[2];
	std::vector<std::vector<XYZ>> dxda(3,std::vector<XYZ>());
	for(int a=0;a<3;++a){
		dxda[a]=std::vector<XYZ>(3,XYZ());
		for(int m=0;m<3;++m){
			dxda[a][m] =dlfda[a].doridx[m] +localx[0]*dlfda[a].daxisdx[0][m]+
					localx[1]*dlfda[a].daxisdx[1][m]+localx[2]*dlfda[a].daxisdx[2][m];
		}
	}
	std::vector<XYZ> dvda(3,XYZ());
	for(int a=0;a<3;++a){
			for(int m=0;m<3;++m){
				for(int k=0;k<3;++k){
					dvda[a][m] += dvdx[k]*dxda[a][m][k];
				}
			}
		}
	return dvda;
}
LocalFrame NSPgeometry::make_localframe(XYZ origin, XYZ px, XYZ pxy ){
	LocalFrame f;
	f.origin_=origin;
	px=px-origin;
	pxy=pxy-origin;
	f.axis_.clear();
	XYZ xa=px/sqrt(px.squarednorm());
	XYZ za=cross(px,pxy);
	za = za/sqrt(za.squarednorm());
	XYZ ya= cross(za,xa);
	f.axis_.push_back(xa);
	f.axis_.push_back(ya);
	f.axis_.push_back(za);
	return f;
}
LocalFrame NSPgeometry::make_localframe(const QuaternionCrd & Q, const XYZ & ori ) {
	LocalFrame lf;
	lf.origin_=ori;
	XYZ axis;
	axis.x_=1.0-2.0*(Q[2]*Q[2]+Q[3]*Q[3]);
	axis.y_=2.0*(Q[1]*Q[2]-Q[0]*Q[3]);
	axis.z_=2.0*(Q[1]*Q[3]+Q[0]*Q[2]);
	lf.axis_.push_back(axis);
	axis.x_=2.0*(Q[2]*Q[1]+Q[0]*Q[3]);
	axis.y_=1.0-2.0*(Q[1]*Q[1]+Q[3]*Q[3]);
	axis.z_=2.0*(Q[2]*Q[3]-Q[0]*Q[1]);
	lf.axis_.push_back(axis);
	axis.x_=2.0*(Q[3]*Q[1]-Q[0]*Q[2]);
	axis.y_=2.0*(Q[3]*Q[2]+Q[0]*Q[1]);
	axis.z_=1.0-2.0*(Q[1]*Q[1]+Q[2]*Q[2]);
	lf.axis_.push_back(axis);
	return lf;
}

LocalFrame NSPgeometry::make_localframe(const Line & zaxis, const XYZ & ponx ){
	LocalFrame f;
	f.origin_=zaxis.perpendicular_foot(ponx);
	XYZ xa=(1/distance(ponx,f.origin_))*(ponx-f.origin_);
	XYZ ya=cross(zaxis.direction,xa);
	f.axis_.clear();
	f.axis_.push_back(xa);
	f.axis_.push_back(ya);
	f.axis_.push_back(zaxis.direction);
	return f;
}
