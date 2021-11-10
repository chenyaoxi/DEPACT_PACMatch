/*
 * testbackboneenergy.cpp
 *
 *  Created on: 2017年8月2日
 *      Author: hyliu
 */
#include "backbone/backboneenergy.h"
#include "backbone/backbonemoves.h"
#include "dstl/randomengine.h"
#include <iostream>
#include <fstream>
using namespace NSPproteinrep;
std::set<int> sspositions(const std::vector<BackBoneSite> &chain) {
	std::set<int> res;
	for (int i = 0; i < chain.size(); ++i) {
		char ss = chain.at(i).sscodechar();
		if (ss == 'H' || ss == 'E')
			res.insert(i);
	}
	return res;
}
int main(int argc, char **argv) {
	std::vector<BackBoneSite> chain;
	readbackbonesites(std::string(argv[1]), chain);
	std::string filename = std::string("movedx.pdb");
	std::ofstream ofs;
	ofs.open(filename.c_str());
	writeSitesToPDB(ofs, chain);
	ofs.close();
	std::set<int> phipsifixed = sspositions(chain);
	std::vector<std::string> controllines;
//	controllines.push_back(std::string("RefPBSeq=FromConf"));
	controllines.push_back(std::string("PBPackingWeight=0.2"));
	initenergycontrols(controllines);
//	initenergycontrols();
	ChainEnergyControl  ce=prepareenergycontrol(&chain,std::string());
	BackBoneEnergy energy;
	double ene = energy.totalenergy(chain, &ce);
	std::cout << "Initial Energy: " << ene << std::endl;
	int seed = std::stoi(std::string(argv[2]));
	NSPdstl::RandomEngine<>::getinstance().reseed(seed);
//	BackBoneMoveSelector ms(chain.size(),phipsifixed);
	BackBoneMoveSelector ms(chain.size());
	BackBoneMoves moves;
	/*	bool cont=true;
	 while(cont){
	 ms.selectmoves(chain,&moves);
	 if(moves.startposi>moves.endposi) {
	 if(moves.nloops()>0)cont=false;
	 }
	 }
	 std::cout<<moves.startposi<<" "<<moves.endposi<<std::endl;
	 energy.deltaE(moves,0);*/
	double temp = 0.5;
	int nstep = 500000;
	for (int step = 0; step < nstep; ++step) {
		ms.selectmoves(chain, &moves);
		int nloops = moves.nloops();
		EnergyComponents pecomp0;
		for (int l = 0; l < nloops; ++l) {
			EnergyComponents pecomp1;
			double de = energy.deltaE(moves, l,&ce, &pecomp0, &pecomp1);
			double ene0 = moves.energy0;
			double paccept = exp(-de / temp);
			if (NSPdstl::RandomEngine<>::getinstance().realrng(0.0, 1.0)()
					< paccept) {
				moves.updatechain(&chain, ce.ssaspbtype(),l);
				moves.energy0 += de;
				ene += de;
				EnergyComponents ecomp=ce.ecomp();
				for (int eidx = 0; eidx < ecomp.energies.size(); ++eidx) {
					ecomp.energies[eidx] += pecomp1.energies[eidx]
							- pecomp0.energies[eidx];
				}
				double detot = 0.0;
				for (int eidx = 1; eidx < ecomp.energies.size(); ++eidx) {
					detot += pecomp1.energies[eidx] - pecomp0.energies[eidx];
				}
				std::cout << "Step: " << de << ": " << detot << std::endl;
				EnergyComponents ecompt;
				ce.ecomp()=ecompt;
				std::cout << "Step: " << step << " Energies: " << ene << " : "
						<< energy.totalenergy(chain, &ce) << std::endl;
				for (int eidx = 0; eidx < ecomp.energies.size(); ++eidx) {
					std::cout << "\t" << ecomp.energies[eidx] << " : "
							<< ce.ecomp().energies[eidx] << std::endl;
				}
				pecomp0.energies = pecomp1.energies;
//				std::cout <<"Step: " <<step<<" Energy diff: "<< ene-energy.totalenergy(chain)<<std::endl;
//				std::cout <<moves.startposi <<" "<< moves.endposi <<std::endl;
//				continue;
			}
		}
		if ((step + 1) % 1 == 0) {
			EnergyComponents ecompsave=ce.ecomp();
			EnergyComponents ecomprecalc;
			ce.ecomp()=ecomprecalc;
			std::cout << "Step: " << step << " Energies: " << ene << " : "
					<< energy.totalenergy(chain, &ce) << std::endl;
			ecomprecalc=ce.ecomp();
			for (int eidx = 0; eidx < ecomprecalc.energies.size(); ++eidx) {
				std::cout << "\t" << ecomprecalc.energies[eidx] << " : "
						<< ecompsave.energies[eidx] << std::endl;
			}
			ce.ecomp()=ecompsave;
			std::string filename = "moved" + std::to_string(step) + ".pdb";
			std::ofstream ofs;
			ofs.open(filename.c_str());
			writeSitesToPDB(ofs, chain);
			ofs.close();
		}
	}
}

