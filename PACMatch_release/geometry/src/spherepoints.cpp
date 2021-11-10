/*
 * genspherepoints.cpp
 *
 *  Created on: 2016年4月2日
 *      Author: hyliu
 */

#include "geometry/spherepoints.h"
#include <cmath>
using namespace NSPgeometry;
void NSPgeometry::genspherepoints( int npoints, std::vector<std::vector<double>> & results,double ra) {
	 double gr = (sqrt(5.0)+1.0)/2.0 - 1.0;
	 double PI=acos(-1.0);
	 double gang=gr*2.0*PI;
	 for (int i=0; i<npoints; i++) {
		 std::vector<double> p(3);
		 p[2]=2.0*(double)(i+1)/(double)npoints-1.0;
		 double r=1.0-p[2]*p[2];
		 if(r<0.0) r=0.0;
		 r=sqrt(r);
		 double phi=gang*(double)(i+1);
		 p[0]=r*cos(phi);
		 p[1]=r*sin(phi);
		 for(auto &c:p) c *=ra;
		 results.push_back(p);
	 }
}
void NSPgeometry::genspherepoints( int npoints, std::vector<XYZ> & results,double ra) {
	 double gr = (sqrt(5.0)+1.0)/2.0 - 1.0;
	 double PI=acos(-1.0);
	 double gang=gr*2.0*PI;
	 for (int i=0; i<npoints; i++) {
		 std::vector<double> p(3);
		 p[2]=2.0*(double)(i+1)/(double)npoints-1.0;
		 double r=1.0-p[2]*p[2];
		 if(r<0.0) r=0.0;
		 r=sqrt(r);
		 double phi=gang*(double)(i+1);
		 p[0]=r*cos(phi);
		 p[1]=r*sin(phi);
		 results.push_back(ra*XYZ(p[0],p[1],p[2]));
	 }
}

void NSPgeometry::genspherepoints(const LocalFrame &lf, int npoints, std::vector<XYZ> & results,double ra) {
	 double gr = (sqrt(5.0)+1.0)/2.0 - 1.0;
	 double PI=acos(-1.0);
	 double gang=gr*2.0*PI;
	 for (int i=0; i<npoints; i++) {
		 std::vector<double> p(3);
		 p[2]=2.0*(double)(i+1)/(double)npoints-1.0;
		 double r=1.0-p[2]*p[2];
		 if(r<0.0) r=0.0;
		 r=sqrt(r);
		 double phi=gang*(double)(i+1);
		 p[0]=r*cos(phi);
		 p[1]=r*sin(phi);
		 results.push_back(lf.local2globalcrd(ra*XYZ(p[0],p[1],p[2])));
	 }
}

