/*
 * comparelooptorsionscores.cpp
 *
 *  Created on: 2017年3月29日
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
	ifs.close();
	ifs.open(argv[3]);   //generated loops
	std::ofstream ofs;
	ofs.open("looptorsionscores.dat");
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
		std::vector<BackBoneSite> chain;
		chain.resize(lp[3]-lp[0],BackBoneSite());
		std::copy(sites.begin()+lp[0],sites.begin()+lp[3],chain.begin());
		std::vector<std::vector<BackBoneSite>> newloops;
		newloops.resize(ncopies+1,std::vector<BackBoneSite>());
		int looplength=lp[2]-lp[1]+2;
		newloops[0].resize(looplength, BackBoneSite());
		std::copy(sites.begin() + lp[1] - 1, sites.begin() + lp[2] + 1,
						newloops[0].begin());
		newloops[0][0].sscode='C';
		newloops[0][looplength-1].sscode='C';
		bool readsucess=true;
		for(unsigned int m=1;m<ncopies+1; ++m) {
			if(!readbackbonesites(ifs,looplength,newloops[m])){
				readsucess=false;
				break;
			}
		}
		if(!readsucess) break;
//		ofs <<"Begin loop "<< nloop<<std::endl;
		if (ncopies<10) continue;
		std::vector<double> scores;
//		if(nloop<=733) continue;
		for (unsigned int i=0; i<ncopies + 1; ++i) {
			std::copy(newloops[i].begin(), newloops[i].begin()+looplength,
					chain.begin()+lp[1]-lp[0]-1);
			BackBoneSite &pns=*(chain.begin()+lp[1]-lp[0]-2);
			BackBoneSite &ns=*(chain.begin()+lp[1]-lp[0]-1);
			ns.phi(pns);
			BackBoneSite &cs=*(chain.begin()+lp[2]-lp[0]);
			BackBoneSite &ncs=*(chain.begin()+lp[2]-lp[0]+1);
			cs.psi(ncs);
//			cs.omiga(ncs);
			double score=0.0;
			for(int i=lp[1]-lp[0]-1; i<= lp[2]-lp[0];++i){
				BackBoneSite &is=chain[i];
				BackBoneSite &isp=chain[i-1];
				BackBoneSite &isn=chain[i+1];
				const PhiPsiDistr *distr=&PhiPsiDistr::coildistr();
				if(is.resname=="GLY") {
					distr=&PhiPsiDistr::glydistr();
				} else if(is.resname =="PRO") {
					if(isp.omiga() >-90.0 && isp.omiga()<90.0) distr=&PhiPsiDistr::cisprodistr();
					else distr=&PhiPsiDistr::transprodistr();
				} else {
					if (isn.resname=="PRO") distr=&PhiPsiDistr::preprodistr();
				}
				score += distr->statisticalenergy(is.phi(),is.psi());
			}
			scores.push_back(score);
		}
		ofs <<nloop;
		for (auto e:scores) {
			ofs <<"\t"<< e;
		}
		ofs <<std::endl;
//		ofs <<"End loop "<< nloop <<std::endl;
	}
}



