/*
 * hbonded.h
 *
 *  Created on: 2016年12月10日
 *      Author: hyliu
 */

#ifndef BACKBONE_HBONDED_H_
#define BACKBONE_HBONDED_H_

#include "backbone/backbonesite.h"

namespace NSPproteinrep{
struct HBondGeometry {
	double rha;   // d-H...a-b , hydrogen-acceptor distance
	double adha;  // donor-hydrogen-acceptor angle
	double ahab;  // hydrogen-acceptor-next_atom angle
	double rda;
	double adab;
};
void segmenthbonded(const std::vector<BackBoneSite> &seg,std::vector<unsigned int> &hbonded);
void updatehbonded(const std::vector<BackBoneSite> &seg1,std::vector<unsigned int> &hbonded1,
		const std::vector<BackBoneSite> &seg2,std::vector<unsigned int> &hbonded2);
unsigned int hbonded(const BackBoneSite &bs1, const BackBoneSite &bs2, HBondGeometry *hbg=nullptr);
bool n1o2hbonded(const BackBoneSite & bs1,const BackBoneSite &bs2,HBondGeometry *hbg=nullptr);
}




#endif /* BACKBONE_HBONDED_H_ */
