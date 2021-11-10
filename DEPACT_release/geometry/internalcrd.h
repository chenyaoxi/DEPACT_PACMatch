/*
 * internalcrd.h
 *
 *  Created on: 2016年11月13日
 *      Author: hyliu
 */

#ifndef GEOMETRY_INTERNALCRD_H_
#define GEOMETRY_INTERNALCRD_H_
#include "geometry/calculators.h"
namespace NSPgeometry{
template <typename ATOMKEY>
struct InternalCrd{
	typedef ATOMKEY AtomKeyType;
	InternalCrd(double d=0.0, double a=0.0, double t=0.0):distance(d),angle(a),torsion(t){;}
	InternalCrd(const AtomKeyType &i, const AtomKeyType &j, const AtomKeyType &k,
			double d=0.0, double a=0.0, double t=0.0):  ia(i),ja(j),ka(k), distance(d),
			angle(a),torsion(t){;}
	void setInternal(double d=0.0, double a=0.0, double t=0.0){
		distance=d;angle=a;torsion=t;}
	template <typename CRDMAP>
	void calcInternal(const CRDMAP & crdmap, const XYZ &rl) {
		if(undefined() ) return;
		distance=NSPgeometry::distance(crdmap[ka],rl);
		angle=NSPgeometry::angle(crdmap[ja],crdmap[ka],rl);
		torsion=NSPgeometry::torsion(crdmap[ia],crdmap[ja],crdmap[ka],rl);
	}
	template <typename CRDMAP>
	XYZ calcXYZ(const CRDMAP & crdmap) const {
		if(undefined()) return XYZ;
		return InternaltoXYZ(crdmap[ka],crdmap[ja],crdmap[ia], distance, angle, torsion);
	}

	bool undefined() const {return ia == ja;}
	double distance;
	double angle;
	double torsion;
	AtomKeyType ia;
	AtomKeyType ja;
	AtomKeyType ka;
};
}


#endif /* GEOMETRY_INTERNALCRD_H_ */
