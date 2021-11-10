/*
 * localframe.h
 *
 *  Created on: 2016年3月31日
 *      Author: hyliu
 */

#ifndef LOCALFRAME_H_
#define LOCALFRAME_H_
#include "geometry/calculators.h"
#include "geometry/line.h"
#include "geometry/quaternioncrd.h"
#include <vector>
namespace NSPgeometry {

/*! Defines a local Cartesian coordinate system using three
 * orthonormal basis vectors starting from an origin point.
 *
 */
struct LocalFrame {

	/*!calculates coordinates in the global coordinate system from local coordinates
	 *
	 */
	XYZ local2globalcrd(const XYZ & p) const {;
		return (origin_ + p.x_*axis_[0] + p.y_*axis_[1] + p.z_*axis_[2]);
	}
	/*!calculates coordinates in the local coordinates system from global coordinates
	 *
	 */
	XYZ global2localcrd(XYZ  p) const {
		p=p-origin_;
		return XYZ(dot(p,axis_[0]), dot(p,axis_[1]),dot(p,axis_[2]));
	}
	/*!a line (a point and a direction) corresponding to axis i
	 *
	 */
	Line axisline(int i) const {return Line(origin_,axis_[i]);}

	/*!three orthonormal basis vectors
	 *
	 */
	std::vector<XYZ> axis_;

	/*!origin of the local coordinate system
	 *
	 */
	XYZ origin_;
};

struct DLocalFrameDx {
	std::vector<XYZ> doridx;
	std::vector<std::vector<XYZ>> daxisdx;
};

struct A3LocalFrame: public LocalFrame {
	A3LocalFrame(){;}
	A3LocalFrame (XYZ a1, XYZ a2,XYZ a3);
	std::vector<XYZ> distributedvdx(const XYZ & x, const XYZ & dvdx) const;
	std::vector<DLocalFrameDx> dlfda;
};
/*!make a local frame with the given origin, a point on the positive local x axis and a point
 * on the first squadron of the local xy plane
 *
 */
LocalFrame make_localframe(XYZ origin, XYZ vec1, XYZ vec2);

/*!make a local frame with a line along the local z axis and a point on the positive local x axis.
 *
 */
LocalFrame make_localframe(const Line &zaxis, const XYZ & pointonx);

/*!make a local frame by rotating the global x, y, z axis with Q and move the origin to ori.
 *
 */
LocalFrame make_localframe(const QuaternionCrd & Q, const XYZ & ori=XYZ(0.0,0.0,0.0) );
}

#endif /* LOCALFRAME_H_ */
