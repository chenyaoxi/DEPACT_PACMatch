/*
 * rigidtransform.cpp
 *
 *  Created on: 2016年11月8日
 *      Author: hyliu
 */

#include "geometry/rigidtransform.h"
#include "dstl/randomengine.h"
#include "geometry/quatfit.h"

using namespace NSPgeometry;

RigidTransform NSPgeometry::applyRotation(const Rotation &r, const RigidTransform &rt){
	RigidTransform rtnew;
	typename Rotation::constmatrixpointer m=rt.rotation().matrix();
	typename Rotation::constmatrixpointer rm=r.matrix();
	Rotation & rnew=rtnew.rotation();
	XYZ & trans=rtnew.translation();
	typename Rotation::matrixpointer mnew=rnew.matrix();
	for(int i=0;i<3;++i) {
		for(int j=0; j<3;++j)
			mnew[i][j]= rm[i][0]*m[0][j]+rm[i][1]*m[1][j]+rm[i][2]*m[2][j];
	}
	rnew.center()= rt.rotation().center();
	trans= rt.translation()+rnew.center();
	r.apply(&trans);
	trans= trans -rnew.center();
	return rtnew;
}

void RigidTransform::apply(XYZ *p) const{
	rot_.apply(p);
	*p = *p + trans_;
}
XYZ RigidTransform::applytoCopy(const XYZ &p) const {
	XYZ res=p;
	apply(&res);
	return res;
}
RigidTransform NSPgeometry::randomrigidtransform(double maxrotate, double maxtranslate){
	static auto &reng=NSPdstl::RandomEngine<>::getinstance();
	double angle=maxrotate*(reng.realrng(0.0,1.0)()-0.5);
	XYZ axis(reng.realrng(0.0,1.0),1.0);
	return RigidTransform(QuaternionCrd(axis,angle),XYZ(0.0,0.0,0.0),
				XYZ(reng.realrng(0.0,1.0),maxtranslate));
}

RigidTransform NSPgeometry::superpose(const std::vector<XYZ> &crda, const std::vector<XYZ> &crdb,
		const std::vector<std::pair<int,int>> & alignedpositions,double *rmsd2){
		std::vector<XYZ> refcrd;
		std::vector<XYZ> crd;
		for(auto & ap:alignedpositions){
			refcrd.push_back(crda[ap.first]);
			crd.push_back(crdb[ap.second]);
		}
		QuatFit qf;
		double dev2=qf.setup(refcrd,crd);
		if(rmsd2) *rmsd2=dev2;
		return qf.getRigidTransform();
}

