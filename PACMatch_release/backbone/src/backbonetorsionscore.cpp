/*
 * backbonetorsionscore.cpp
 *
 *  Created on: 2016年12月15日
 *      Author: hyliu
 */

#include "backbone/backbonetorsionscore.h"
#include "pdbstatistics/phipsidistr.h"

using namespace NSPproteinrep;

double NSPproteinrep::backbonetorsionscore(const std::vector<BackBoneSite> &loop) {
	bool iscis=false;
	double score=0.0;
	for(unsigned int i=0;i<loop.size();++i) {
		std::string resname=loop[i].resname;
		const NSPpdbstatistics::PhiPsiDistr *distr=&(NSPpdbstatistics::PhiPsiDistr::coildistr());
		if(resname=="GLY") {
			distr=&(NSPpdbstatistics::PhiPsiDistr::glydistr());
		} else if(resname=="PRO") {
			if(iscis) distr=&(NSPpdbstatistics::PhiPsiDistr::cisprodistr());
			else distr=&(NSPpdbstatistics::PhiPsiDistr::transprodistr());
			iscis=false;
			if(loop[i].omiga() >-90.0 && loop[i].omiga() <90.0) iscis=true;
		}
		if(i<loop.size()-1 && resname != "GLY" && resname != "PRO") {
			if(loop[i+1].resname == "PRO") {
				distr=&(NSPpdbstatistics::PhiPsiDistr::preprodistr());
			}
		}
		score +=distr->statisticalenergy(loop[i].phi(),loop[i].psi());
	}
	return score;
}

