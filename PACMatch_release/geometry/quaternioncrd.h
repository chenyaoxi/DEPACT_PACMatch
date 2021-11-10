/*
 * quaterioncrd.h
 *
 *  Created on: 2016年3月30日
 *      Author: hyliu
 */

#ifndef QUATERNIONCRD_H_
#define QUATERNIONCRD_H_
#include "calculators.h"
#include <iostream>

#include <vector>
#include <cmath>
namespace NSPgeometry {

struct LocalFrame;
/*!A quaternion represents a rotation around the origin.
 *
 */
struct QuaternionCrd{
	/*!default constructor generates an identical rotation
	 *
	 */
	QuaternionCrd(): Q({1.0,0.0,0.0,0.0}){;}

	/*!A random rotation in uniformly distributed  quaternion space.
	 * rng() should return a uniform random number between 0 to 1.0
	 * dumy is an arbitrary integer to resolve function overloading
	 */
	template <typename RNG>
	QuaternionCrd(RNG &rng,int dumy){
	bool get=false;
	while (!get) {
		double norm=0.0;
		std::vector<double> q;
		for(int d=0; d<4; ++d) {
			double x=rng();  //uniform random number between 0 to 1.0
			x=x*2.0-1.0;
			q.push_back(x);
			norm +=x*x;
		}
		if(norm >1.0) continue;
		get=true;
		norm=1.0/sqrt(norm);
		for(int d=0; d<4; ++d) Q.push_back(q[d]*norm);
		}
	}

	/*!The rotation from the global coordinate system to local coordinate system lf
	 *
	 */
	QuaternionCrd(const LocalFrame & lf);

	static double scale_rotation (double angle) {
		if(angle <0) angle =-angle;
		double u=(angle - sin(angle))/3.14159265;
		return pow(u,1.0/3.0);
	}

	QuaternionCrd(const std::vector<double> & q) {
		Q=q;
	}

	/*!rotation around xyz by angle (in degree by default)
	 *
	 */
	QuaternionCrd(XYZ xyz,double angle,double deg=3.14159265358979323846/180.0) {
		xyz=xyz/sqrt(xyz.squarednorm());
		Q.resize(4);
		angle= deg*angle;
		Q[0]=cos(angle*0.5);
		double sintheta=sin(angle*0.5);
		if(Q[0] < 0) {
			Q[0]=-Q[0];
			sintheta=-sintheta;
		}
		Q[1]=xyz.x_*sintheta;
		Q[2]=xyz.y_*sintheta;
		Q[3]=xyz.z_*sintheta;
	}

	/*!get rotation angle in degree
	 *
	 */
	double angle() {
//		std::cout <<Q[0]<<"  ";
		double a=2*acos(Q[0]);
		a=a*180.0/3.14159265;
		if(a>180.0) a -=360.0;
		return a;
	}

	/*!map a rotation (point on the surface of a 4-d sphere) to a point
	 * within a 3d-sphere
	 *
	 */
	void maptoR3(XYZ &s, XYZ &simage) {
		s=XYZ(Q[1],Q[2],Q[3]);
		s=s/sqrt(s.squarednorm());
		double a=angle_rad();
		if(a<0) {s= -1*s; a=-a;}
		s=scale_rotation(a)*s;
		a=2*3.14159265-a;
		simage=scale_rotation(a)*s;
	}
	/*!get rotation angle in radius
	 *
	 */
	double angle_rad() {
//		std::cout <<Q[0]<<"  ";
		double a=2*acos(Q[0]);
		if(a>3.14159265) a -=2*3.14159265;
		return a;
	}
	XYZ axis() {
		XYZ axis(Q[1],Q[2],Q[3]);
		double n=sqrt(axis.squarednorm());
		if(n>0.0000000001) return axis/n;
		return XYZ(0.0,0.0,0.0);
	}
	XYZ vectcomp() const {
		return XYZ(Q[1],Q[2],Q[3]);
	}
	/*!get the rotation matrix
	 *
	 */
	void getmatrix(double (&matrix)[3][3]) const;

	/*!difference in the 4-d space
	 *
	 */
	static std::vector<double> diff(std::vector<double> q1, std::vector<double> q2) {
		double n1=0.0,n2=0.0;
		std::vector<double> tmp1,tmp2;
		for(int d=0; d<4;++d) {
			tmp1[d]=q1[d]-q2[d];
			n1 += tmp1[d]*tmp1[d];
			tmp2[d]=q1[d]+q2[d];
			n2 += tmp2[d]*tmp2[d];
		}
		if(n1< n2) return tmp1;
		return tmp2;
	}

	/*!q1 and -q1 represents the same rotation.
	 * one of q1 and -q1, which is closer to q2 in the 4-d space, is returned.
	 *
	 */
	static std::vector<double> shift(std::vector<double> q1, const std::vector<double> q2) {
		double n1=0.0,n2=0.0;
		std::vector<double> qequ(4);
		for(int d=0; d<4; ++d) qequ[d]=-q1[d];
		double tmp;
		for(int d=0; d<4;++d) {
			tmp=q1[d]-q2[d];
			n1 += tmp*tmp;
			tmp=qequ[d]-q2[d];
			n2 += tmp*tmp;
		}
		if(n1< n2) return q1;
		return qequ;
	}

	static std::vector<double> mean(std::vector<double> sum, double  count) {
		double n=0.0;
		for(int d=0; d<4;++d) {
			sum[d] = sum[d]/count;
			n+= sum[d]*sum[d];
		}
		n=sqrt(n);
		for (int d=0; d<4; ++d) sum[d] = sum[d]/n;
		return sum;
	}

	/*!
	 * return reverse rotation (the same angle, inverted axis)
	 */
	QuaternionCrd invert() {
		QuaternionCrd res;
		res[0]=Q[0];
		for (int i=1; i<4;++i)	res[i]=-Q[i];
		return res;
	}

	/*!
	 * -q is quivalent to q
	 */
	QuaternionCrd negative_equivalent() {
		QuaternionCrd res;
		for(int i=0; i<4; ++i) res[i]=-Q[i];
		return res;
	}

	const double & operator[] (int i) const {
		return Q[i];
	}
	double &operator[] (int i) { return Q[i];}
	/*!
	 * The angle of the rotation changing Q2 to this, in degree
	 */
	double diff(QuaternionCrd Q2);
	/*!
	 * The angle of the rotation changing Q2 to this, in radius
	 */
	double diff_rad(QuaternionCrd Q2);
	/*!
	 * The 4 components of the quaternion
	 */
	std::vector<double> Q;
};
/*!
 * combined rotation Q1Q2
 */
inline QuaternionCrd operator *(const QuaternionCrd & q1, const QuaternionCrd &q2){
	QuaternionCrd Q;
	XYZ v1=q1.vectcomp();
	XYZ v2=q2.vectcomp();
	Q[0]=q1[0]*q2[0]-dot(v1,v2);
	XYZ v=q1[0]*v2 + q2[0]*v1 + cross(v1,v2);
	Q[1]=v.x_;
	Q[2]=v.y_;
	Q[3]=v.z_;
	return Q;
}
/**
 * calculate the Quaternion that rotate the axis "from" to the axis "to
 */
QuaternionCrd quaternionaligntwovectors(const XYZ & to, const XYZ & from);
}



#endif /* QUATERNIONCRD_H_ */
