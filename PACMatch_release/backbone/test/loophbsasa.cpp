/*
 * loophbsasa.cpp
 *
 *  Created on: 2016年12月12日
 *      Author: hyliu
 */

#include "backbone/backbonesasa.h"
#include "backbone/segments.h"
#include "backbone/hbonded.h"
#include <iostream>
#include <fstream>
using namespace NSPproteinrep;

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
	for(auto & lp: loops) {
		std::vector<std::vector<BackBoneSite>> segments;
		std::vector<BackBoneSite> loop;
		segments.push_back(std::vector<BackBoneSite>());
		segments[0].resize(lp[1] - lp[0] - 1, BackBoneSite());
		std::copy(sites.begin() + lp[0], sites.begin() + lp[1] - 1,
				segments[0].begin());
		segments.push_back(std::vector<BackBoneSite>());
		segments[1].resize(lp[3] - lp[2] - 1, BackBoneSite());
		std::copy(sites.begin() + lp[2] + 1, sites.begin() + lp[3],
				segments[1].begin());
		loop.resize(lp[2] - lp[1] + 2, BackBoneSite());
		std::copy(sites.begin() + lp[1] - 1, sites.begin() + lp[2] + 1,
				loop.begin());

		std::vector<std::vector<BackBoneSASA>> segmentsasa;
		std::vector<BackBoneSASA> loopsasa;
		for (auto & seg : segments) {
			segmentsasa.push_back(std::vector<BackBoneSASA>());
			segmentSASA(seg, segmentsasa.back());
		}
		segmentSASA(loop, loopsasa);
		updatesegmentSASA(segmentsasa[0], segmentsasa[1]);
		updatesegmentSASA(segmentsasa[0], loopsasa);
		updatesegmentSASA(segmentsasa[1], loopsasa);

		std::vector<unsigned int> loophbonded;
		loophbonded.resize(loop.size(), 0u);
		unsigned int nhb = 0;
		for (unsigned int m = 0; m < loop.size(); ++m) {
			for (unsigned int i = 0; i < 2; ++i) {
				for (unsigned int n = 0; n < segments[i].size(); ++n) {
					unsigned int hb = hbonded(loop[m], segments[i][n]);
					if (hb)
						++nhb;
					if (hb == 1) {
						loophbonded[m] += 1;
					} else if (hb == 2) {
						loophbonded[m] += 10;
					}
				}
			}
			if (m == loop.size() - 1)
				continue;
			for (unsigned int n = m + 1; n < loop.size(); ++n) {
				unsigned int hb = hbonded(loop[m], loop[n]);
				if (hb)
					++nhb;
				if (hb == 1) {
					loophbonded[m] += 1;
					loophbonded[n] += 10;
				} else if (hb == 2) {
					loophbonded[m] += 10;
					loophbonded[n] += 1;
				}
			}
		}

		for(unsigned int i=0;i<loop.size();++i) {
			std::vector<double>exposed;
			loopsasa[i].getexposed(&exposed);
			unsigned int nhbonded=loophbonded[i]%10;
			unsigned int ohbonded=loophbonded[i]/10u;
			std::cout <<loop[i].resname<<"N: " << nhbonded <<" " <<exposed[0] <<std::endl;
			std::cout <<loop[i].resname<<"O: " << ohbonded <<" " <<exposed [3]  <<std::endl;
		}
	}
}
