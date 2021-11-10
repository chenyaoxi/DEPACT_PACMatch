/*
 * structalign.h
 *
 *  Created on: 2017年8月20日
 *      Author: hyliu
 */

#ifndef GEOMETRY_STRUCTALIGN_H_
#define GEOMETRY_STRUCTALIGN_H_
#include "geometry/calculators.h"
#include "dstl/alignset.h"
#include <vector>
namespace NSPgeometry{

inline NSPdstl::AlignedPositions alignpointset(const std::vector<XYZ> &seta, const std::vector<XYZ> &setb,
		double rcut){
	double rcut2=rcut*rcut;
	auto matcher=[rcut2](const XYZ &p1, const XYZ &p2)->bool{
		return (p1-p2).squarednorm()<rcut2;
	};
	return NSPdstl::SetMatch::alignset(seta,setb,matcher);
}
}



#endif /* GEOMETRY_STRUCTALIGN_H_ */
