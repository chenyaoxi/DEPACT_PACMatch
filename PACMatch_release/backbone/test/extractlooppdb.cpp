/*
 * extractlooppdb.cpp
 *
 *  Created on: 2017年1月3日
 *      Author: hyliu
 */

#include "backbone/backboneloop.h"
#include "backbone/backbonesasa.h"
#include "backbone/hbonded.h"
#include "backbone/torsionvectorscorer.h"
#include <iostream>
#include <fstream>
using namespace NSPproteinrep;

int main(int argc, char **argv) {
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
//	std::vector<BackBoneSite> tmpsites;
//	readbackbonesites("tmplatesites.dat",tmpsites);
//	TorsionVectorScorer &tvscorer=TorsionVectorScorer::getinstance(&tmpsites);
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
	unsigned int loopid=std::stoi(std::string(argv[4]));
	std::ofstream ofs;
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
/*		std::vector<std::vector<BackBoneSite>> segments;
		segments.push_back(std::vector<BackBoneSite>());
		segments[0].resize(lp[1] - lp[0] - 1, BackBoneSite());
		std::copy(sites.begin() + lp[0], sites.begin() + lp[1] - 1,
				segments[0].begin());
		segments.push_back(std::vector<BackBoneSite>());
		segments[1].resize(lp[3] - lp[2] - 1, BackBoneSite());
		std::copy(sites.begin() + lp[2] + 1, sites.begin() + lp[3],
				segments[1].begin());*/
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
			for(auto it=newloops[m].begin(); it!=newloops[m].end();++it){
				it->chainid=newloops[0][0].chainid;
			}
		}
		if(!readsucess) break;
		if(nloop!= loopid) continue;
		for (unsigned int i=0; i<ncopies + 1; ++i) {
			std::copy(newloops[i].begin(), newloops[i].begin()+looplength,
					chain.begin()+lp[1]-lp[0]-1);
			std::string ofile=chain[0].pdbid+"loop"+std::to_string(lp[1]-lp[0])+"_"
				+std::to_string(i)+".pdb";
			ofs.open(ofile.c_str());
			writeSitesToPDB(ofs,chain);
			ofs.close();
		}
		break;
	}
}



