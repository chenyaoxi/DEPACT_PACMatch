/*
 * backbonedecoys.cpp
 *
 *  Created on: 2017年8月6日
 *      Author: hyliu
 */
#include "backbone/backboneenergy.h"
#include "backbone/backbonemoves.h"
#include "dstl/randomengine.h"
#include <iostream>
#include <fstream>
#include <sstream>
using namespace NSPproteinrep;
/*void writedecoy(std::ostream &os, int decoyid,
		const std::vector<BackBoneSite> &chain) {
	os << decoyid << " " << chain[0].pdbid << " " << chain.size() << std::endl;
	for (auto s : chain) {
		if (s.sscode == ' ')
			s.sscode = 'x';
		if (s.resname == "CISPRO")
			s.resname = "PRO";
		os << s.toString();
	}
}
bool readnextdecoy(std::istream &is, int *decoyid,
		std::vector<BackBoneSite> *chain) {
	std::string pdbid;
	int chainsize;
	char buffer[120];
	int nsection = 0;
	is.getline(buffer, 120);
	if (!is.good())
		return false;
	std::string line = std::string(buffer);
	std::stringstream sstr(line);
	sstr >> *decoyid >> pdbid >> chainsize;
//	if(!sstr.eof() && !sstr.good()) return false;
	chain->clear();
	return readbackbonesites(is, chainsize, *chain);
}
std::set<int> sspositions(const std::vector<BackBoneSite> &chain) {
	std::set<int> res;
	for (int i = 0; i < chain.size(); ++i) {
		char ss = chain.at(i).sscodechar();
		if (ss == 'H' || ss == 'E')
			res.insert(i);
	}
	return res;
}*/
int generatedecoys(std::vector<BackBoneSite> &chain,int nsolutions,int ntries,std::ostream & os) {
	std::vector<LoopMover> movers=getloopmovers(chain,5,15);
	if(movers.empty()) return 0.0;
	ChainEnergyControl ce=prepareenergycontrol(&chain,std::string());
	BackBoneEnergy energy;
	BackBoneMoves move;
	int decoyid=0;
	for(auto &mover:movers) {
		mover.makemoves(chain,nsolutions,ntries, &move);
		int nloops=move.nloops();
		EnergyComponents ecomp0;
		for(int i=0;i<nloops;++i){
			EnergyComponents ecomp1;
			energy.deltaE(move,i,&ce,&ecomp0,&ecomp1);
			std::vector<double> de(5,0.0);
			for(int i=0;i<5;++i) de[i]=ecomp1.energies[i]-ecomp0.energies[i];
			if(de[EnergyComponents::CLASH]>0) continue;
//			if(de[EnergyComponents::BLOCKLOCAL] > 0.0) continue;
			os <<de[EnergyComponents::PHIPSI] <<" "<<de[EnergyComponents::BLOCKLOCAL]<<
														" "<<de[EnergyComponents::BLOCKPACKING]<<std::endl;
			++decoyid;
			if(NSPdstl::RandomEngine<>::getinstance().realrng(0.0,1.0)()< 0.01){
				std::vector<BackBoneSite> nchain=chain;
				move.updatechain(&nchain,ce.ssaspbtype(),i);
				std::string pdb=chain[0].pdbid+".pdb";
				std::string pdb1=chain[0].pdbid+"_decoy"+std::to_string(decoyid)+".pdb";
				std::ofstream ofs1;
				ofs1.open(pdb.c_str());
				writeSitesToPDB(ofs1,chain);
				ofs1.close();
				ofs1.open(pdb1.c_str());
				writeSitesToPDB(ofs1,nchain);
				EnergyComponents ecompa,ecompb;
				ce.ecomp()=ecompa;
				energy.totalenergy(chain,&ce);
				ecompa=ce.ecomp();
				ce.ecomp()=ecompb;
				energy.totalenergy(nchain,&ce);
				ecompb=ce.ecomp();
				std::vector<double> de2(5,0.0);
				for(int i=0;i<5;++i) de2[i]=ecompb.energies[i]-ecompa.energies[i];
				for(int i=0;i<5;++i) {
					std::cout <<" "<<de[i] <<":"<<de2[i];
				}
				std::cout <<std::endl;
			}
		}
	}
}
/*	std::vector<std::string> scheme0 { "PhiPsiWeight=1.0", "PBLocalWeight=0.0",
			"ClashWeight=0.0", "PBPackingWeight=0.0" };
	std::vector<std::string> scheme1 { "PhiPsiWeight=0.0", "PBLocalWeight=0.0",
			"ClashWeight=1.0", "PBPackingWeight=0.0" };
	std::set<int> phipsifixed = sspositions(chain);
	BackBoneMoveSelector ms(chain.size(), phipsifixed);
	int nflexible = chain.size() - phipsifixed.size();
	if (nflexible < 7)
		return 0;
	double phipsitol = 0.005 * (double) nflexible;
	initenergycontrols(scheme0);

	double initphipsi = energy0.totalenergy(chain);
	std::cout << "Initial PhiPsi Energy: " << initphipsi << std::endl;
	adjustenergycontrols(scheme1);
	BackBoneEnergy energy1;
	double temp = 1.0;
	for (int ngenerated = 0; ngenerated < ndecoys; ++ngenerated) {
		std::vector<BackBoneSite> nchain = chain;
		BackBoneMoves moves;
		int movesteps = 0;
		bool done = false;
		while (!done) {
			ms.selectmoves(nchain, &moves);
			int nloops = moves.nloops();
			int incmove = 0;
			for (int l = 0; l < nloops; ++l) {
				double de = energy1.deltaE(moves, l);
				double paccept = exp(-de / temp);
				if (NSPdstl::RandomEngine<>::getinstance().realrng(0.0, 1.0)()
						< paccept) {
					std::vector<BackBoneSite> trialchain = nchain;
					moves.updatechain(&trialchain, l);
					double deltaphipsi = energy0.totalenergy(trialchain)
							- initphipsi;
					if (deltaphipsi < phipsitol) {
//						std::cout << "Delta PhiPsi Energy: " << deltaphipsi
//								<< std::endl;
						moves.updatechain(&nchain, l);
						moves.energy0 += de;
						incmove = 1;
						break;
					}
				}
			}
			movesteps += incmove;
			if (movesteps >= nmoves) {
				double phipsi = energy0.totalenergy(nchain);
				if (phipsi <= initphipsi + phipsitol) {
					writedecoy(os, decoyid++, nchain);
					done = true;
				}
			}
		}
	}
	return ndecoys;
}*/

int main(int argc, char**argv) {
	std::string filename = std::string("deltaE.dat");
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
	if(argc > 2) {
		unsigned int seed=std::stoi(std::string(argv[2]));
		NSPdstl::RandomEngine<>::getinstance().reseed(seed);
	}
	std::ofstream ofs;
	ofs.open(filename.c_str());
	int maxchains = 500;
	int nchains = 0;
	auto itstart = sites.begin();
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
//				std::string pdbfile = chain[0].pdbid + "_ref" + ".pdb";
//				std::ofstream ofspdb;
//				ofspdb.open(pdbfile.c_str());
//				writeSitesToPDB(ofspdb, chain);
//				ofspdb.close();
				generatedecoys(chain, 500, 500, ofs);
				++nchains;
			}
		}
		itstart = itend + 1;
	}
	ofs.close();
/*	std::ifstream ifs;
	ifs.open(filename.c_str());
	for (int i = 0; i < ndecoys; ++i) {
		std::vector<BackBoneSite> chain;
		int decoyid;
		readnextdecoy(ifs, &decoyid, &chain);
		std::string pdbfile = chain[0].pdbid + "_decoy_"
				+ std::to_string(decoyid) + ".pdb";
		std::ofstream ofs;
		ofs.open(pdbfile.c_str());
		writeSitesToPDB(ofs, chain);
		ofs.close();
	}*/
}
