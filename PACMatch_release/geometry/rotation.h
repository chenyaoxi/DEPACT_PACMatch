/*
 * rotation.h
 *
 *  Created on: 2015年11月9日
 *      Author: hyliu
 */

#ifndef ROTATION_H_
#define ROTATION_H_
#include <cassert>
#include "geometry/xyz.h"
#include "geometry/quaternioncrd.h"

namespace NSPgeometry {

/*! Defines a rotation around an axis passing a point.
 *
 * For a position represented by a column vector xyz, the rotation operation
 * leads to this transformed position: matrix_*(xyz - center_) + center_
 */
class Rotation {

public:
	Rotation(){;}

	/*!Rotation around a 3d space vector from p2 to p1 by phi.
	 *
	 * Default p2 is the origin. p1-p2 distance must not be zero.
	 */
//	void init(double phi, const XYZ & p1, const XYZ &p2 = XYZ());

	/*!The Rotation that when applied to rl, the torsion angle ri-rj-rk-rl becomes phi
	 *
	 */
	void init(const XYZ & ri, const XYZ &rj, const XYZ & rk, const XYZ &rl,
			double phi);

	Rotation(const XYZ & ri, const XYZ &rj, const XYZ & rk, const XYZ &rl,
			double phi) {
		init(ri,rj,rk,rl,phi);
	}

	void init(const QuaternionCrd & Q, const XYZ &center){
		Q.getmatrix(matrix_);
		center_=center;
	}
	Rotation(const QuaternionCrd & Q, const XYZ &center){
		init(Q,center);
	}
	/*!Rotate and change the input point
	 *
	 */
	void apply(XYZ *rin) const {
		XYZ t = *rin - center_;
		rin->x_ = matrix_[0][0] * t.x_ + matrix_[0][1] * t.y_
				+ matrix_[0][2] * t.z_ + center_.x_;
		rin->y_ = matrix_[1][0] * t.x_ + matrix_[1][1] * t.y_
				+ matrix_[1][2] * t.z_ + center_.y_;
		rin->z_ = matrix_[2][0] * t.x_ + matrix_[2][1] * t.y_
				+ matrix_[2][2] * t.z_ + center_.z_;
	}

	/*!Rotate and returns the position. Input position unchanged
	 *
	 */
	XYZ applytoCopy(const XYZ & rin) const {
		XYZ t = rin;
		apply(&t);
		return t;
	}

	XYZ &center() {return center_;}
	const XYZ &center() const {return center_;}
	typedef double (*matrixpointer) [3];
	typedef const double (*constmatrixpointer)[3];
	matrixpointer matrix() {return matrix_;}
	constmatrixpointer matrix() const {return matrix_;}

private:
	/*! rotation matrix
	 */
	double matrix_[3][3] {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};

	/*! A point on the axis
	 *
	 */
	XYZ center_ {0.0,0.0,0.0};
};
Rotation rotationaligntwovectors(const XYZ &tolocation,const XYZ &todirection, const XYZ &fromdirection);

}

#endif /* ROTATION_H_ */
