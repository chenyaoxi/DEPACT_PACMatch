/*
 * nativevsdecoyloops.cpp
 *
 *  Created on: 2017年10月12日
 *      Author: hyliu
 */
#include "backbone/backboneenergy.h"
#include "backbone/closealoop.h"
#include "dstl/randomengine.h"
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
	unsigned int nloop=0;

	std::vector<std::string> controllines;
//	controllines.push_back(std::string("RefPBSeq=FromConf"));
	initenergycontrols(controllines);
	int lpid=0;
	for (auto &lp:loops) {
		if(lpid<=306) {
			lpid++;
			continue;
		}
		std::vector<BackBoneSite> chain;
		std::vector<BackBoneSite> nativeloop;
		chain.resize(lp[3]-lp[0]);
		nativeloop.resize(lp[2]-lp[1]);
		std::copy(sites.begin()+lp[0],sites.begin()+lp[3],chain.begin());
		std::copy(sites.begin()+lp[1],sites.begin()+lp[2],nativeloop.begin());
		int loopstart=lp[1]-lp[0];
		int looplength=lp[2]-lp[1];
		ChainEnergyControl  ce=prepareenergycontrol(&chain,std::string());
		BackBoneEnergy energy;
		energy.partialE(chain,loopstart,loopstart+looplength,&ce);
		EnergyComponents enative=ce.ecomp();
		std::cout <<lpid<<" native-energies:";
		for(auto e:enative.energies) std::cout <<" "<< e;
		double ephipsi_native=enative.energies[EnergyComponents::PHIPSI];
		int ndecoys=0;
		int ntry=0;
		bool done=false;
		while (!done) {
			CloseALoop clp(chain,loopstart,looplength);
			auto decoys=clp.getsolutions();
			for(auto decoy:decoys){
				std::copy(decoy->begin(),decoy->begin()+looplength, chain.begin()+loopstart);
				ChainEnergyControl  ced=prepareenergycontrol(&chain,std::string());
				energy.partialE(chain,loopstart,loopstart+looplength,&ced);
				EnergyComponents edecoy=ced.ecomp();
				if(edecoy.energies[EnergyComponents::TOTAL] >=9000) continue;
				double dephipsi=edecoy.energies[EnergyComponents::PHIPSI]-ephipsi_native;
				if(dephipsi >= 0.1*(double) looplength) continue;
				ndecoys++;
				std::cout << " decoy-energies:";
				for(auto e:edecoy.energies) std::cout <<" " <<e;
				if(ndecoys>=20) {
					std::cout <<" done";
					done=true;
					break;
				}
			}
			ntry++;
			if(ntry>=10000) done=true;
		}
		std::cout<<std::endl;
		lpid++;
	}
}
