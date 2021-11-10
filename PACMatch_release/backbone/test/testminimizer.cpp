/*
 * testminimizer.cpp
 *
 *  Created on: 2017年8月11日
 *      Author: hyliu
 */

#include "backbone/backboneminimizer.h"
#include "dstl/randomengine.h"

using namespace NSPproteinrep;
using namespace NSPdstl;

int main(int argc, char**argv) {

	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
	if (argc > 2) {
		unsigned int seed = std::stoi(std::string(argv[2]));
		NSPdstl::RandomEngine<>::getinstance().reseed(seed);
	}
	int maxchains = 10;
	int nchains = 0;
	auto itstart = sites.begin();
	MinimizerControl mct;
	mct.maxsteps=500000;
	mct.maxnochangesteps=10000;
	mct.settriphase(1.0,0.2,std::vector<int>({1001,1501,5001}));
	mct.nprint=100;
	while (nchains < maxchains && itstart < sites.end()) {
		auto itend = itstart;
		while (!chainendsite(itend, sites.end()))
			++itend;
		int chainlength = itend - itstart + 1;
		if (chainlength >= 150 && chainlength < 500) {
			if (fragstartsite(itstart, sites.end(), chainlength, std::string(),
					false)) {
				std::vector<BackBoneSite> chain(chainlength);
				std::copy(itstart, itstart + chainlength, chain.begin());
				std::string filename=chain[0].pdbid+"_init.pdb";
				std::ofstream ofs;
				ofs.open(filename.c_str());
				writeSitesToPDB(ofs,chain);
				BackBoneMinimizer min;
				std::shared_ptr<std::vector<BackBoneSite>> minchain=min.run(mct,chain);
				++nchains;
			}
		}
		itstart = itend + 1;
	}
}

