/*
 * gengeratedvsnativeloops.cpp
 *
 *  Created on: 2016年12月14日
 *      Author: hyliu
 */

#include "backbone/backboneloop.h"
#include "backbone/backbonesasa.h"
#include "backbone/hbonded.h"
#include "backbone/backbonetorsionscore.h"
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
	ifs.close();
	ifs.open(argv[3]);   //generated loops
	std::ofstream ofs;
	ofs.open("loopproperties.dat");
	while (ifs.good()) {
		char buffer[120];
		ifs.getline(buffer,120);
		std::string line(buffer);
		if( !ifs.good()) break;
		std::stringstream sstr(line);
		unsigned int nloop, ncopies;
		sstr>> nloop >>ncopies;
		if(ncopies == 0) continue;
		std::vector<long> &lp=loops[nloop];
		std::vector<std::vector<BackBoneSite>> segments;
		segments.push_back(std::vector<BackBoneSite>());
		segments[0].resize(lp[1] - lp[0] - 1, BackBoneSite());
		std::copy(sites.begin() + lp[0], sites.begin() + lp[1] - 1,
				segments[0].begin());
		segments.push_back(std::vector<BackBoneSite>());
		segments[1].resize(lp[3] - lp[2] - 1, BackBoneSite());
		std::copy(sites.begin() + lp[2] + 1, sites.begin() + lp[3],
				segments[1].begin());
		std::vector<std::vector<BackBoneSASA>> segmentsasa;
		std::vector<std::vector<unsigned int>> segmenthbonded;
		for (auto & seg : segments) {
			segmentsasa.push_back(std::vector<BackBoneSASA>());
			segmentSASA(seg, segmentsasa.back());
			segmenthbonded.push_back(std::vector<unsigned int>());
			NSPproteinrep::segmenthbonded(seg,segmenthbonded.back());
		}
		updatesegmentSASA(segmentsasa[0], segmentsasa[1]);
		updatehbonded(segments[0],segmenthbonded[0],segments[1],segmenthbonded[1]);

		std::vector<std::vector<BackBoneSite>> newloops;
		newloops.resize(ncopies+1,std::vector<BackBoneSite>());
		int looplength=lp[2]-lp[1]+2;
		newloops[0].resize(looplength, BackBoneSite());
		std::copy(sites.begin() + lp[1] - 1, sites.begin() + lp[2] + 1,
						newloops[0].begin());
		bool readsucess=true;
		for(unsigned int m=1;m<ncopies+1; ++m) {
			if(!readbackbonesites(ifs,looplength,newloops[m])){
				readsucess=false;
				break;
			}
		}
		if(!readsucess) break;
		ofs <<"Begin loop "<< nloop<<std::endl;
		for (unsigned int i=0; i<ncopies + 1; ++i) {
			ofs <<i<<"Torsion_score: " << NSPproteinrep::backbonetorsionscore(newloops[i])<<std::endl;
			std::vector<std::vector<BackBoneSASA>> sasa_f=segmentsasa;
			std::vector<std::vector<unsigned int>> hbonded_f=segmenthbonded;
			std::vector<BackBoneSASA> loopsasa;
			std::vector<unsigned int> loophbonded;
			segmentSASA(newloops[i],loopsasa);
			NSPproteinrep::segmenthbonded(newloops[i],loophbonded);
			updatesegmentSASA(sasa_f[0], loopsasa);
			updatesegmentSASA(sasa_f[1], loopsasa);
			updatehbonded(segments[0],hbonded_f[0],newloops[i],loophbonded);
			updatehbonded(segments[1],hbonded_f[1],newloops[i],loophbonded);
			for(unsigned int s=0; s<2; ++s){
				for(unsigned int p=0; p<segments[s].size();++p) {
					std::vector<double> exposed;
					std::vector<double> exposed_c;
					sasa_f[s][p].getexposed(&exposed);
					segmentsasa[s][p].getexposed(&exposed_c);
					if(exposed[0] != exposed_c[0]) {
						ofs <<i<<" "<<"N_expose_change: " <<segmenthbonded[s][p]%10 <<" " <<hbonded_f[s][p]%10 <<" "
								<<exposed_c[0] <<" "<< exposed[0] <<std::endl;
					}
					if(exposed[3] != exposed_c[3]) {
						ofs <<i<<" "<<"C_expose_change: " <<segmenthbonded[s][p]/10u <<" " <<hbonded_f[s][p]/10u <<" "
								<<exposed_c[3] <<" "<< exposed[3] <<std::endl;
					}
				}
			}
			for(unsigned int p=0;p<newloops[i].size();++p) {
					std::vector<double>exposed;
					loopsasa[p].getexposed(&exposed);
					unsigned int nhbonded=loophbonded[p]%10;
					unsigned int ohbonded=loophbonded[p]/10u;
					if(newloops[i][p].resname != "PRO")
						ofs <<i<<" "<<"N: " << nhbonded <<" " <<exposed[0] <<std::endl;
					ofs <<i<<" "<<"O: " << ohbonded <<" " <<exposed [3]  <<std::endl;
			}
		}
		ofs <<"End loop "<< nloop <<std::endl;
	}
}
