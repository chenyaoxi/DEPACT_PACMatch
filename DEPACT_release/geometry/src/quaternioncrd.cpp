/*
 * quaternioncrd.cpp
 *
 *  Created on: 2016年4月1日
 *      Author: hyliu
 */

#ifndef SRC_QUATERNIONCRD_CPP_
#define SRC_QUATERNIONCRD_CPP_

#include "geometry/quaternioncrd.h"
#include "geometry/localframe.h"
#include <cmath>

using namespace NSPgeometry;
double QuaternionCrd::diff(QuaternionCrd Q2) {
	Q2=(*this)*Q2.invert();
	return Q2.angle();
}
double QuaternionCrd::diff_rad(QuaternionCrd Q2) {
	Q2=(*this)*Q2.invert();
	return Q2.angle_rad();
}

void QuaternionCrd::getmatrix(double (&R)[3][3]) const {
		R[0][0]=1.0-2.0*(Q[2]*Q[2]+Q[3]*Q[3]);
		R[0][1]=2.0*(Q[1]*Q[2]-Q[0]*Q[3]);
		R[0][2]=2.0*(Q[1]*Q[3]+Q[0]*Q[2]);

		R[1][0]=2.0*(Q[2]*Q[1]+Q[0]*Q[3]);
		R[1][1]=1.0-2.0*(Q[3]*Q[3]+Q[1]*Q[1]);
		R[1][2]=2.0*(Q[2]*Q[3]-Q[0]*Q[1]);

		R[2][0]=2.0*(Q[3]*Q[1]-Q[0]*Q[2]);
		R[2][1]=2.0*(Q[3]*Q[2]+Q[0]*Q[1]);
		R[2][2]=1.0-2.0*(Q[1]*Q[1]+Q[2]*Q[2]);
}
QuaternionCrd::QuaternionCrd(const LocalFrame & lf) {
	  std::vector<std::vector<double>> a;
	  for (int i=0; i<3; i++) {
		  std::vector<double> b(3);
		  b[0]=lf.axis_[i].x_;
		  b[1]=lf.axis_[i].y_;
		  b[2]=lf.axis_[i].z_;
		  a.push_back(b);
	  }
	  double trace = a[0][0] + a[1][1] + a[2][2];
	  Q.resize(4);
	  if( trace > 0 ) {
	    double s = 0.5 / sqrt(trace+ 1.0);
	    Q[0] = 0.25 / s;
	    Q[1] = ( a[2][1] - a[1][2] ) * s;
	    Q[2] = ( a[0][2] - a[2][0] ) * s;
	    Q[3] = ( a[1][0] - a[0][1] ) * s;
	  } else {
	    if ( a[0][0] > a[1][1] && a[0][0] > a[2][2] ) {
	      double s = 2.0 * sqrt( 1.0 + a[0][0] - a[1][1] - a[2][2]);
	      Q[0] = (a[2][1] - a[1][2] ) / s;
	      Q[1] = 0.25 * s;
	      Q[2] = (a[0][1] + a[1][0] ) / s;
	      Q[3] = (a[0][2] + a[2][0] ) / s;
	    } else if (a[1][1] > a[2][2]) {
	      double s = 2.0 * sqrt( 1.0 + a[1][1] - a[0][0] - a[2][2]);
	      Q[0] = (a[0][2] - a[2][0] ) / s;
	      Q[1] = (a[0][1] + a[1][0] ) / s;
	      Q[2] = 0.25f * s;
	      Q[3] = (a[1][2] + a[2][1] ) / s;
	    } else {
	      double s = 2.0 * sqrt( 1.0f + a[2][2] - a[0][0] - a[1][1] );
	      Q[0] = (a[1][0] - a[0][1] ) / s;
	      Q[1] = (a[0][2] + a[2][0] ) / s;
	      Q[2] = (a[1][2] + a[2][1] ) / s;
	      Q[3] = 0.25 * s;
	    }
	  }
}

QuaternionCrd NSPgeometry::quaternionaligntwovectors(const XYZ & to, const XYZ & from){
	static double epsinon=1.e-16;
	double rto=to.squarednorm();
	double rfrom=from.squarednorm();
	if(rto <= epsinon || rfrom <=epsinon) return QuaternionCrd();
	double angle=acos(dot(to,from)/sqrt(rto*rfrom))*180.0/3.14159265;
	XYZ axis=cross(from,to);
	if(axis.squarednorm() <= 1.e-16) return QuaternionCrd();
	return QuaternionCrd(axis,angle);
}

#endif /* SRC_QUATERNIONCRD_CPP_ */
