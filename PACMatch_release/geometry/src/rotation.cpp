/*
 * rotation.cpp
 *
 *  Created on: 2015年11月9日
 *      Author: hyliu
 */
#include "geometry/calculators.h"
#include "geometry/rotation.h"
#include <cmath>
#include <exception>
using namespace NSPgeometry;
#ifndef PI
#define PI 3.1415926535897932
#endif

void Rotation::init(const XYZ & ri, const XYZ &rj, const XYZ & rk, const XYZ &rl,double phi) {
	double p=phi - torsion(ri, rj, rk, rl);
	init(QuaternionCrd(rk-rj,p*180.0/PI),rk);
}

Rotation NSPgeometry::rotationaligntwovectors(const XYZ &tolocation,const XYZ &todirection,  const XYZ &fromdirection){
	Rotation R;
	QuaternionCrd q=quaternionaligntwovectors(todirection,fromdirection);
	R.init(q,tolocation);
	return R;
}

