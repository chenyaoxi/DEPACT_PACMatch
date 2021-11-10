/*
 * backbonesasa.h
 *
 *  Created on: 2016年12月10日
 *      Author: hyliu
 */

#ifndef BACKBONE_BACKBONESASA_H_
#define BACKBONE_BACKBONESASA_H_
#include "geometry/atomsasa.h"
#include "backbonesite.h"

namespace NSPproteinrep {

class BackBoneSASA {
public:
	static std::vector<NSPgeometry::SASAParameter> backboneSASAparameters;
	enum {N,CA,C,O,CB,CD_PRO};
	BackBoneSASA(){;}
	BackBoneSASA(const BackBoneSite &bs);
	void init(const BackBoneSite &bs);
	void reset(const BackBoneSite &bs);
	void getexposed(std::vector<double> *exposed){
		for(auto &a:atomsasa_)exposed->push_back(a.exposedfraction());
	}
	friend void updateSASA(BackBoneSASA &s1, BackBoneSASA &s2);
private:
	std::vector<NSPgeometry::AtomSASA> atomsasa_;
};
void updateSASA(BackBoneSASA &s1, BackBoneSASA &s2);

void segmentSASA(const std::vector<BackBoneSite> &segment,std::vector<BackBoneSASA> &sasa);
void updatesegmentSASA(std::vector<BackBoneSASA> &sasa1, std::vector<BackBoneSASA> &sasa2);
}


#endif /* BACKBONE_BACKBONESASA_H_ */
