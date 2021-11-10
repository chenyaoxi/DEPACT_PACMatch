/*
 * xyz.h
 *
 *  Created on: 2015年11月9日
 *      Author: hyliu
 */

#ifndef XYZ_H_
#define XYZ_H_
#include "stdio.h"
#include <iostream>
#include <string>
#include <vector>
namespace NSPgeometry {
/*! A vector or point in Cartesian space
 *
 */
struct XYZ {
	double x_;
	double y_;
	double z_;
	/*!initializer, default initialized to origin
	 *
	 */
	XYZ(double ix = 0.0, double iy = 0.0, double iz = 0.0) :
			x_(ix), y_(iy), z_(iz) {
		;
	}

	/*! p[0],p[1],and p[2] mapped to x_,y_ and z_, respectively
	 *
	 */
	XYZ (const std::vector<double> & p) : x_(p[0]), y_(p[1]),z_(p[2]){;}
	double squarednorm() const {
		return x_ * x_ + y_ * y_ + z_ * z_;
	}

	/*! Draw a random point from a uniform distribution
	 * within a sphere of radius ra around the origin.
	 *  rng() should return a (random) number between 0 and 1.
	 */
	template <typename RNG>
	XYZ (RNG &rng,double ra){
		bool get=false;
		while (!get) {
			double norm=0.0;
			std::vector<double> p;
			for(int d=0; d<3; ++d) {
				double x=rng();  //uniform random number between 0 to 1.0
				x=x*2.0-1.0;
				p.push_back(x);
				norm +=x*x;
			}
			if(norm >1.0) continue;
			get=true;
			x_=p[0]*ra;
			y_=p[1]*ra;
			z_=p[2]*ra;
		}
	}
	XYZ operator-() const {
		return XYZ(-x_,-y_,-z_);
	}

	/*!maps a vector of XYZ objects to a (three-times-longer) vector of doubles
	 * The new doubles are APPENDED to the crd vector.
	 */
	static void xyzstovector(const std::vector<XYZ> & points, std::vector<double> & crd) {
		for(const XYZ & p:points) {
			crd.push_back(p.x_);
			crd.push_back(p.y_);
			crd.push_back(p.z_);
		}
	}
	/*!maps a vector of doubles to a (three-times-shorted) vector of XYZ objects
	 * the new XYZ objects are APPENDED to the res vector.
	 */
	static void vectortoxyzs(const std::vector<double> &crd, std::vector<XYZ> &res){
		for(int i=0; i<crd.size(); i+=3){
			res.push_back(XYZ(crd[i],crd[i+1],crd[i+2]));
		}
	}
	/*!maps a XYZ object to a vector of three doubles
	 *
	 */
	std::vector<double> tovector() const {
		std::vector<double> p(3);
		p[0]=x_; p[1]=y_;p[2]=z_;
		return p;
	}
	double & operator[] (int i) { if(i==0) return x_; else if(i==1) return y_; else if(i==2) return z_;
	else {std::cout <<"XYZ index out of range" <<std::endl;abort();}}
	const double & operator[] (int i) const { if(i==0) return x_; else if(i==1) return y_; else if(i==2) return z_;
	else {std::cout <<"XYZ index out of range" <<std::endl;abort();}}

	/*!
	 * map to a formatted string
	 */
	std::string toString() const {
		const std::string FMT {"%8.3f%8.3f%8.3f"};
		char s[30];
		sprintf(s, FMT.c_str(), x_, y_, z_);
		return (std::string(s));
	}

	/*!
	 * check if the coordinates are within an appropriate range.
	 */
	bool valid(double min=-1.0e10, double max=1.0e10) const {
		return (x_ > min && x_ <max) &&
			   (y_>min && y_<max) &&
			   (z_>min && z_<max);
	}
};
inline XYZ operator -(const XYZ &p1, const XYZ &p2) {
	XYZ p;
	p.x_ = p1.x_ - p2.x_;
	p.y_ = p1.y_ - p2.y_;
	p.z_ = p1.z_ - p2.z_;
	return p;
}
inline XYZ operator +(const XYZ &p1, const XYZ &p2) {
	XYZ p;
	p.x_ = p1.x_ + p2.x_;
	p.y_ = p1.y_ + p2.y_;
	p.z_ = p1.z_ + p2.z_;
	return p;
}

inline XYZ operator *(const XYZ &p, double d) {
	XYZ r;
	r.x_ = p.x_ * d;
	r.y_ = p.y_ * d;
	r.z_ = p.z_ * d;
	return r;
}
inline XYZ operator *(double d, const XYZ &p) {
	return p * d;
}
inline XYZ operator /(const XYZ &p, double d) {
	return p * (1.0 / d);
}
inline double dot(const XYZ &p1, const XYZ &p2) {
	return p1.x_ * p2.x_ + p1.y_ * p2.y_ + p1.z_ * p2.z_;
}
inline XYZ cross(const XYZ &p1, const XYZ &p2) {
	XYZ res;
	res.x_ = p1.y_ * p2.z_ - p1.z_ * p2.y_;
	res.y_ = p1.z_ * p2.x_ - p1.x_ * p2.z_;
	res.z_ = p1.x_ * p2.y_ - p1.y_ * p2.x_;
	return res;
}

} //namespace geometry

#endif /* XYZ_H_ */
