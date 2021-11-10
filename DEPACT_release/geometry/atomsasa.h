/*
 * atomsasa.h
 *
 *  Created on: 2016年12月10日
 *      Author: hyliu
 */

#ifndef GEOMETRY_ATOMSASA_H_
#define GEOMETRY_ATOMSASA_H_
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include <set>
#include <vector>
namespace NSPgeometry {

struct SASAParameter {
	double solventradius;
	double soluteradius;
	double nsurfacepoints;
	SASAParameter(double ru=1.4,double rv=1.4, int nsurf=256): soluteradius(ru),solventradius(rv),nsurfacepoints(nsurf){;}
};

class AtomSASA {
public:
	AtomSASA(){;}
	AtomSASA( const SASAParameter &par):par_(par){;}
	AtomSASA(const LocalFrame &lf, const SASAParameter &par);
	void init(const LocalFrame &lf, const SASAParameter &par);
	void init(const LocalFrame &lf);
	double exposedfraction() {return (double)(exposed_.size())/(double) (surfacepoints_.size());}
	XYZ &crd() {return solutecrd_;}
	const XYZ &crd() const {return solutecrd_;}
	const SASAParameter & par() const {return par_;}
	void update(const XYZ &crdb, double rb);
private:
	XYZ solutecrd_;
	std::vector<XYZ> surfacepoints_;
	std::set<unsigned int> exposed_;
	SASAParameter par_;
};

void updateAtomSASA (AtomSASA & atoma, AtomSASA &atomb);
std::vector<double> calc_sasa(const std::vector<XYZ> &crd,
		const std::vector<double> &radii_solu, double radius_solv,int npoint=256);

}



#endif /* GEOMETRY_ATOMSASA_H_ */
