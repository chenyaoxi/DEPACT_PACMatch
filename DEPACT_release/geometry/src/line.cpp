/*
 * line.cpp
 *
 *  Created on: 2016年11月3日
 *      Author: hyliu
 */

#include "geometry/line.h"

using namespace NSPgeometry;

XYZ Line::perpendicular_foot(const XYZ & point) const{
	return origin+dot(point-origin,direction)*direction;
}
bool Line::contains(const XYZ &point) const {
	static const double EPSINON=1.e-20;
	if(distancefrom (point) <= EPSINON){
		return true;
	}
	return false;
}
Line NSPgeometry::make_line(const XYZ & p1, const XYZ &p2){
	return Line(p1,(1/NSPgeometry::distance(p1,p1))*(p2-p1));
}

