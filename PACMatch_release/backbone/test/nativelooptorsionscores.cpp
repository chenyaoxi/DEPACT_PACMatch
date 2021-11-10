/*
 * nativelooptorsionscores.cpp
 *
 *  Created on: 2017年3月30日
 *      Author: hyliu
 */

#include "backbone/backboneloop.h"
#include "pdbstatistics/phipsidistr.h"
#include <iostream>
#include <fstream>
using namespace NSPproteinrep;
using namespace NSPpdbstatistics;

int main(int argc, char **argv) {
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
	std::vector<std::vector<long>> loops;
	std::ifstream ifs;
	ifs.open(argv[2]); // sites files
	while (ifs.good()) {
		long l0, l1, l2, l3;
		ifs >> l0 >> l1 >> l2 >> l3;
		if (ifs.good()) {
			loops.push_back(std::vector<long>( { l0, l1, l2, l3 }));
		}
	}
	std::ofstream ofs;
	ofs.open("looplocations_ts.dat");
	for(int nloop=0; nloop<loops.size();++nloop) {
		std::vector<long> &lp = loops[nloop];
		std::vector<BackBoneSite> chain;
		chain.resize(lp[3] - lp[0], BackBoneSite());
		std::copy(sites.begin() + lp[0], sites.begin() + lp[3], chain.begin());
		int looplength = lp[2] - lp[1] + 2;
		double score = 0.0;
		for (int i = lp[1] - lp[0] - 1; i <= lp[2] - lp[0]; ++i) {
			BackBoneSite &is = chain[i];
			BackBoneSite &isp = chain[i - 1];
			BackBoneSite &isn = chain[i + 1];
			const PhiPsiDistr *distr = &PhiPsiDistr::coildistr();
			if (is.resname == "GLY") {
				distr = &PhiPsiDistr::glydistr();
			} else if (is.resname == "PRO") {
				if (isp.omiga() > -90.0 && isp.omiga() < 90.0)
					distr = &PhiPsiDistr::cisprodistr();
				else
					distr = &PhiPsiDistr::transprodistr();
			} else {
				if (isn.resname == "PRO")
					distr = &PhiPsiDistr::preprodistr();
			}
			score += distr->statisticalenergy(is.phi(), is.psi());
		}
		ofs << lp[0] << '\t' << lp[1] << '\t' << lp[2] << '\t' << lp[3] << '\t'
				<< score << std::endl;
	}
}

