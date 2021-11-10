/*
 * calculators.h
 *
 *  Created on: 2016年2月12日
 *      Author: hyliu
 */

#ifndef CALCULATORS_H_
#define CALCULATORS_H_
#include "geometry/xyz.h"
#include <cmath>

#define NULLANGLE -10000
namespace NSPgeometry {

inline double distance(const XYZ & p1, const XYZ & p2) {
	return std::sqrt((p1 - p2).squarednorm());
}
/*!calculate distance and derivatives of distance over coordinates
 *
 */
double distance(const XYZ & p1, const XYZ & p2,std::vector<XYZ> *deriv);
inline double distance2(const XYZ &p1, const XYZ &p2) {
	return (p1-p2).squarednorm();
}
double angle(const XYZ &p1, const XYZ &p2, const XYZ & p3);

/*!calculate cosine of the angle and derivatives of cosine with respect to coordinates
 *
 */
double cos_angle(const XYZ &p1, const XYZ &p2, const XYZ & p3,std::vector<XYZ> *deriv);
double torsion(const XYZ &p1, const XYZ &p2, const XYZ &p3, const XYZ &p4);

/*!calculate torsional angle (in radius) and its derivatives with respect coordinates
 *
 */
double torsion(const XYZ &p1, const XYZ &p2, const XYZ &p3, const XYZ &p4,std::vector<XYZ> *deriv);
double rmsd(const std::vector<XYZ> & crd1, const std::vector<XYZ> & crd2);
template<typename T> XYZ make_XYZ(const T & t) {
	return XYZ(t.x, t.y, t.z);
}
template<typename T> double distance(const T & t1, const T & t2) {
	XYZ p1 = make_XYZ(t1);
	XYZ p2 = make_XYZ(t2);
	return distance(p1, p2);
}
template<typename T> double angle(const T & t1, const T & t2, const T & t3) {
	XYZ p1 = make_XYZ(t1);
	XYZ p2 = make_XYZ(t2);
	XYZ p3 = make_XYZ(t3);
	return angle(p1, p2, p3);
}
template<typename T> double torsion(const T & t1, const T & t2, const T & t3,
		const T & t4) {
	XYZ p1 = make_XYZ(t1);
	XYZ p2 = make_XYZ(t2);
	XYZ p3 = make_XYZ(t3);
	XYZ p4 = make_XYZ(t4);
	return torsion(p1, p2, p3, p4);
}

/*!calculate rmsd squared between two set of coordinates contained in two containers.
 * the first set are in the range [aiter, aend).
 * The second set are in the range [biter,biter+(aend-aiter)).
 * w contains weight for each XYZ coordinate.
 */
template <typename ITER>
double meansquareddist2 (ITER aiter, ITER aend, ITER biter,
		std::vector<double> w=std::vector<double>()) {
		if(w.empty()) for(int i=0; i<aend-aiter; ++i) w.push_back(1.0);
		double dis2_tot=0.0;
		double wtot=0.0;
		auto witer=w.begin();
		while (aiter != aend) {
			dis2_tot += (*witer)*distance2(*aiter,*biter);
			wtot += *(witer++);
			++aiter;
			++biter;
		}
		return dis2_tot/wtot;
}
/*! generate Cartesian coordinate of atom l, with given bond length, angle and torsion
 *   from three atoms
 *
 */
XYZ InternaltoXYZ(const XYZ & rk, const XYZ &rj, const XYZ & ri, double b,
		double theta, double phi);

/*! generate Cartesian coordinate of atom l, with given bond length and angle to previous two atoms
 * The plane containing rj,rk and rl will be arbitrary.
 */
XYZ InternaltoXYZ(const XYZ & rk, const XYZ &rj, double b, double theta);

/*! generate a Cartesian point rl from rk with distance b,
 * The vector from rk to rl arbitrarily set along the x direction
 */
XYZ InternaltoXYZ(const XYZ & rk, double b);
XYZ center(const std::vector<XYZ> & points,std::vector<double> weights=std::vector<double>());
double radiusgyr(const std::vector<NSPgeometry::XYZ> & crd);
double radiusgyr(const std::vector<double> & crd);
}

#endif /* CALCULATORS_H_ */
