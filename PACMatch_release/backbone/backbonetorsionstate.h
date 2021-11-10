/*
 * backbonetorsionstate.h
 *
 *  Created on: 2016年12月27日
 *      Author: hyliu
 */

#ifndef BACKBONETORSIONSTATE_H_
#define BACKBONETORSIONSTATE_H_
#include "backbone/backbonesite.h"
namespace NSPproteinrep {
class BackBoneTorsionState{

public:
	enum {H,E,C1,C2,C3,PROC1,PROC2,PROC3,GLYC1,GLYC2,GLYC3};
	static int torsionstate(double phi, double psi, std::string residuetype=""){
		int state=C1; //alpha region
		phi=shifttorsion(phi);
		psi=shifttorsion(psi);
		if(phi>0 && (residuetype!="PRO" && residuetype!="Pro" && residuetype!="pro"))state=C3;
		else {
			if(phi>-40.0) {
				if(psi >0) state=C2;
			} else {
				if(phi > (-1.5*psi-40.0)) state=C2;
				if(phi <(-1.5*psi-220.0)) state=C2;
			}
		}
		if(residuetype=="PRO" || residuetype=="Pro" || residuetype=="pro") return 3+state;
		if(residuetype=="GLY"||residuetype=="Gly" ||residuetype=="gly") return 6+state;
		return state;
	}
	static int torsionstate(const BackBoneSite &bs){
		if(bs.sscode=='H') return H;
		if(bs.sscode=='E') return E;
		return torsionstate(bs.phi(),bs.psi(),bs.resname);
	}
private:
	static double shifttorsion(double t) {
		while(t<-180.0) t+=360.0;
		while(t>180.0) t-=360.0;
		return t;
	}
};

}
#endif /* BACKBONETORSIONSTATE_H_ */
