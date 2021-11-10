/*
 * chainpackminimizer.cpp
 *
 *  Created on: 2017年8月17日
 *      Author: hyliu
 */


#include "backbone/chainpackminimizer.h"
#include "backbone/backboneenergy.h"
#include "backbone/backbonemoves.h"
#include "dstl/randomengine.h"

using namespace NSPproteinrep;
using namespace NSPpdbstatistics;
using namespace NSPdstl;
void ChainPackMinimizer::writeminconf(const std::string &filename){
	std::ofstream ofs;
	std::string sites=filename+".sites";
	ofs.open(sites.c_str());
	int resseq=1;
	std::vector<BackBoneSite> concat;
	for(auto &chain:*minconf_){
		for(auto &s:chain) {
			if(s.isgap) continue;
			if (s.sscode == ' ')
				s.sscode = 'x';
			if (s.resname == "CISPRO")
				s.resname = "PRO";
			s.resseq=resseq;
			s.resid=resseq++;
			ofs << s.toString();
			concat.push_back(s);
		}
	}
	ofs.close();
	std::string pdb=filename +".pdb";
	ofs.open(pdb.c_str());
	ofs<< "Emin: " <<emin_<<std::endl;
	writeSitesToPDB(ofs, concat);
	ofs.close();
}
std::shared_ptr<std::vector<std::vector<BackBoneSite>>> ChainPackMinimizer::run(
		const MinimizerControl & control, const std::vector<std::vector<BackBoneSite>> &chain) {
//	Chain localchain = chain;
	ChainPack  cpk(chain);
//	std::vector<std::string> energycontrolines { { "RefPBSeq=" }, {
//			"PBLocalWeight=0.0" }, { "PBPackingWeight=0.0" } };
/*	std::vector<std::string> energycontrolines {
		{"PBLocalWeight=1.0" },
		{ "PBPackingWeight=0.3" },
	{"PositionMask=MaskLoopPositions"},
	{"RefPBSeq=FromConf"}
	};
*/
	std::vector<std::string> energycontrolines {
		{"PBLocalWeight=1.0" },
		{ "PBPackingWeight=0.3" },
	{"RefPBSeq=FromConf"},
	{"SSasPBType=1"}
	};
	adjustenergycontrols(energycontrolines);
	cpk.init();
	BackBoneEnergy energy;
	std::string optname("chainpack_opt");
	EnergyComponents ecomp;
	double ene = energy.totalenergy(cpk, &ecomp);
	e0_ = ene;
	CpckMoveSelector ms;
	ChainPackMoves moves;
	auto & rng=NSPdstl::RandomEngine<>::getinstance();
	double temp=0.0;
	bool stop=false;
	while (!stop) {
		if (nstep_ % control.nprint == 0) {
			EnergyComponents ecomprecalc;
			double enerecalc = energy.totalenergy(cpk, &ecomprecalc);
			std::cout << "Step: " << nstep_ <<"Temp: "<< temp << " Energies: ";
			std::cout <<ene <<":" <<enerecalc;
			int tidx=0;
			for (auto e : ecomp.energies)
				std::cout << "\t" <<ecomprecalc.energies[tidx++];
			std::cout << std::endl;
			if(minconf_) writeminconf(optname);
		}
		ms.makemoves(cpk, &moves);
		int nsol=moves.nmoves();
		temp = control.gettemperature(nstep_);
		for (int l = 0; l < nsol; ++l) {
			double de = energy.deltaE(cpk,&moves,l);
			double paccept=exp(-de / temp);
			if (rng.realrng(0.0,1.0)()< paccept) {
				moves.updatechainpack(&cpk, l);
				ene += de;
			}
/*			EnergyComponents ecomprecalc;
			double enerecalc = energy.totalenergy(cpk, &ecomprecalc);
			if(ene-enerecalc >1.e-3 || ene-enerecalc<-1.e-3) {
				std::cout <<"something wrong" <<std::endl;
			}
*/
			if (ene < emin_) {
				emin_ = ene;
				nochangesteps_ = 0;
				minconf_ = std::shared_ptr <std::vector<std::vector<BackBoneSite>> >
				(new std::vector<std::vector<BackBoneSite>>(cpk));
			}
		}
		++nstep_;
		++nochangesteps_;
		stop = control.stop(nstep_, nochangesteps_);
	}
	if (minconf_) {
		std::cout << "Minimization done. E0: " << e0_ << " Emin: " << emin_
				<< std::endl;
		writeminconf(optname);
	} else {
		std::cout << "Minimization failed, no lower energy configuration found."
				<< std::endl;
	}
	return minconf_;
}


