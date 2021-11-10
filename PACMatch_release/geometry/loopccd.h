/*
 * loopccd.h
 *
 *  Created on: 2016年11月23日
 *      Author: hyliu
 */

#ifndef GEOMETRY_LOOPCCD_H_
#define GEOMETRY_LOOPCCD_H_
#include "geometry/line.h"
namespace NSPgeometry{
double ccdangle(const Line &axis, const std::vector<XYZ> &fixpoints, const std::vector<XYZ> &mvpoints);
double rotationtormsd(const Line &axis,double rotation, const std::vector<XYZ> &fixpoints, std::vector<XYZ> mvpoints);
}


#endif /* GEOMETRY_LOOPCCD_H_ */
