/*
 * atomsasa.cpp
 *
 *  Created on: 2016年12月10日
 *      Author: hyliu
 */

#include "geometry/atomsasa.h"
#include "geometry/calculators.h"
#include "geometry/spherepoints.h"

using namespace NSPgeometry;

AtomSASA::AtomSASA(const LocalFrame &lf, const SASAParameter &par){
	init(lf,par);
}
void AtomSASA::init(const LocalFrame &lf, const SASAParameter &par){
	par_=par;
	init(lf);
}
void AtomSASA::init(const LocalFrame &lf){
	solutecrd_=lf.origin_;
	surfacepoints_.clear();
	double r=par_.soluteradius+par_.solventradius;
	genspherepoints(lf,par_.nsurfacepoints,surfacepoints_,r);
	exposed_.clear();
	for(unsigned int i=0;i<surfacepoints_.size();++i){
		exposed_.insert(i);
	}
}
void AtomSASA::update(const XYZ &crdb, double rb){
	double rb2=(rb+par_.solventradius);
	rb2 =rb2*rb2;
	std::vector<unsigned int> toerase;
	for(auto p:exposed_){
		double r2=(surfacepoints_[p]-crdb).squarednorm();
		if(r2<rb2) toerase.push_back(p);
	}
	for(auto e:toerase) exposed_.erase(e);
}
void NSPgeometry::updateAtomSASA(AtomSASA &atoma, AtomSASA &atomb) {
	XYZ ra=atoma.crd();
	XYZ rb=atomb.crd();
	double rab2=(ra-rb).squarednorm();
	double rcut2=atoma.par().soluteradius+2.0*atoma.par().solventradius+atomb.par().soluteradius;
	if(rab2 >rcut2*rcut2) return;
	atoma.update(atomb.crd(),atomb.par().soluteradius);
	atomb.update(atoma.crd(),atoma.par().soluteradius);
}
std::vector<double> NSPgeometry::calc_sasa(const std::vector<XYZ> &crd,
		const std::vector<double> &radii_solu, double radius_solv,int npoints){
	int natoms=crd.size();
	std::vector<double> result(natoms);
	std::vector<AtomSASA> sasa(natoms);
	for(int i=0;i<natoms;++i){
		LocalFrame lf=make_localframe(crd[i],crd[i]+XYZ(1,0,0),crd[i]+XYZ(0,1,0));
		SASAParameter par(radii_solu[i],radius_solv,npoints);
		sasa[i].init(lf,par);
	}
	for(int i=0;i<natoms-1;++i){
		for(int j=i+1;j<natoms; ++j){
			updateAtomSASA(sasa[i],sasa[j]);
		}
	}
	double pi4=4.0*3.14159265;
	for(int i=0;i<natoms;++i){
		double ra=radii_solu[i]+radius_solv;
		result[i]=sasa[i].exposedfraction()*pi4*ra*ra;
	}
	return result;
}


