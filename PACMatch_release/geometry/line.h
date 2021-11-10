/*
 * line.h
 *
 *  Created on: 2016年11月3日
 *      Author: hyliu
 */

#ifndef LINE_H_
#define LINE_H_
#include "geometry/calculators.h"
#include <cmath>
//#include "geometry/localframe.h"
namespace NSPgeometry {
struct Line {
	XYZ origin;
	XYZ direction;
	Line(const XYZ & ori, const XYZ &dir): origin(ori),direction(dir){;}
	XYZ perpendicular_foot(const XYZ &outpoint) const;
	double distancefrom(const XYZ &point) const {return sqrt((point-perpendicular_foot(point)).squarednorm());}
	bool contains(const XYZ &point) const;
};
Line make_line(const XYZ & p1, const XYZ &p2);
}



#endif /* LINE_H_ */
