/*
 * spherepoints.h
 *
 *  Created on: 2016年4月2日
 *      Author: hyliu
 */

#ifndef SPHEREPOINTS_H_
#define SPHEREPOINTS_H_
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include <vector>
namespace NSPgeometry{
void genspherepoints(int npoints, std::vector<std::vector<double>> & results,double r=1.0);
void genspherepoints(int npoints, std::vector<XYZ> & results, double r=1.0);
void genspherepoints(const LocalFrame &lf, int npoints,std::vector<XYZ> &results,double r=1.0);
}

#endif /* SPHEREPOINTS_H_ */
