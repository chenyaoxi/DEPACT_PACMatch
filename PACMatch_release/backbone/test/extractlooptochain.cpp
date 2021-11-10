/*
 * testenergyfunctions.cpp
 *
 *  Created on: 2017年4月25日
 *      Author: hyliu
 */


#include "backbone/backboneloop.h"
#include "backbone/mainchain.h"
#include "backbone/energyfunctions.h"
#include <iostream>
#include <fstream>
#include <algorithm>
using namespace NSPproteinrep;

int main(int argc, char **argv) {
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
//	EnergyTerms eterms;
//	eterms.addtorsionene();
//	eterms.addtetrasefene();
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
//	std::ofstream ofs;
//	ofs.open("loopsefene.dat");
	int loopid=std::stoi(std::string(argv[4]));
	while (ifs.good()) {
		char buffer[120];
		ifs.getline(buffer,120);
		std::string line(buffer);
		if( !ifs.good()) break;
		std::stringstream sstr(line);
		unsigned int nloop, ncopies;
		sstr>> nloop >>ncopies;
		if(nloop != loopid) continue;
		if(ncopies == 0) continue;
		std::vector<long> &lp=loops[nloop];
		MainChain chain;
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
		int loopstart=lp[1]-lp[0]-1;

		for (unsigned int i=0; i<ncopies + 1; ++i) {
			std::copy(newloops[i].begin(),newloops[i].end(),chain.begin()+loopstart);
			chain.setrigidbetween(0,loopstart);
			chain.setrigidbetween(loopstart+looplength-1, chain.size());
			std::string chainname="chain"+std::to_string(loopid)+"_"+std::to_string(i)+".dat";
			std::ofstream ofs;
			ofs.open(chainname.c_str());
			chain.write(ofs);
			ofs.close();
			chainname="chain"+std::to_string(loopid)+"_"+std::to_string(i)+".pdb";
			ofs.open(chainname.c_str());
			writeSitesToPDB(ofs,chain);
		}
		break;
	}
}


