/*
 * backboneminimizer.cpp
 *
 *  Created on: 2017年8月11日
 *      Author: hyliu
 */
#include "backbone/backboneminimizer.h"
#include "backbone/backboneenergy.h"
#include "backbone/backbonemoves.h"
#include "dstl/randomengine.h"

using namespace NSPproteinrep;
using namespace NSPpdbstatistics;
using namespace NSPdstl;
void BackBoneMinimizer::writeminconf(const std::string &filename){
	std::ofstream ofs;
	std::string pdb=filename +".pdb";
	ofs.open(pdb.c_str());
	ofs<< "Emin: " <<emin_<<std::endl;
	writeSitesToPDB(ofs, *minconf_);
	ofs.close();
	std::string sites=filename+".sites";
	ofs.open(sites.c_str());
	for(auto s:*minconf_){
		if(s.isgap) continue;
		if (s.sscode == ' ')
			s.sscode = 'x';
		if (s.resname == "CISPRO")
			s.resname = "PRO";
		ofs << s.toString();
	}
	ofs.close();
}

std::shared_ptr<BackBoneMinimizer::Chain> BackBoneMinimizer::run(
		const MinimizerControl & control, const Chain &chain) {
//	Chain localchain = chain;
	Chain localchain=loops2gaps(chain,1);
//	std::vector<std::string> energycontrolines { { "RefPBSeq=" }, {
//			"PBLocalWeight=0.0" }, { "PBPackingWeight=0.0" } };
/*	std::vector<std::string> energycontrolines {
		{"PBLocalWeight=1.0" },
		{ "PBPackingWeight=0.3" },
	{"PositionMask=MaskLoopPositions"},
	{"RefPBSeq=FromConf"}
	};*/
	std::vector<std::string> energycontrolines {
		{"PBLocalWeight=1.0" },
		{ "PBPackingWeight=0.3" },
	{"RefPBSeq=FromConf"}
	};
	std::string parasetname1("paraset1");
	adjustenergycontrols(parasetname1,energycontrolines);
	ChainEnergyControl ce1=prepareenergycontrol(&localchain,parasetname1);
	std::vector<std::string> energycontrolines2 {
		{"PBLocalWeight=1.0" },
		{ "PBPackingWeight=1.0" },
	{"RefPBSeq=FromConf"}
	};
	std::string set2("paraset2");
	adjustenergycontrols(set2,energycontrolines2);
	ChainEnergyControl ce=prepareenergycontrol(&localchain,set2);
	BackBoneEnergy energy;
	std::string pdbid = chain[0].pdbid;
	std::string optname = pdbid + "_opt";
	double ene = energy.totalenergy(localchain, &ce);
	e0_ = ene;
	bool stop = false;
	BackBoneMoveSelector ms1=make_perturbloopselector(localchain,BackBoneMoveSelector::NONE);
//	BackBoneMoveSelector ms1=make_segmentmoveselector(localchain,BackBoneMoveSelector::NONE);
//	ms1.maxflexposi=7;
//	ms1.maxrotate=50.0;
	BackBoneMoveSelector ms2=make_segmentmoveselector(localchain,BackBoneMoveSelector::ALLNONGAP);
//	BackBoneMoveSelector ms2=make_perturbloopselector(localchain,BackBoneMoveSelector::SSELEMENTS);
	ms2.maxrotate=15.0;
	ms2.maxtranslate=1.5;
	BackBoneMoves moves;
	BackBoneMoveSelector *ms;
	auto & rng=NSPdstl::RandomEngine<>::getinstance();
	double temp=0.0;
	while (!stop) {
		if (nstep_ % control.nprint == 0) {
			std::cout << "Step: " << nstep_ <<"Temp: "<< temp << " Energies: ";
			for (auto e : ce.ecomp().energies)
				std::cout << "\t" << e;
			std::cout << std::endl;
			if(minconf_) writeminconf(optname);
		}
		int nloops = 0;
		int ntry = 0;
		if( rng.realrng(0.0,1.0)()<0.5) ms=&ms1;
		else ms=&ms2;
		while (nloops == 0 && ntry++ < 1000) {
			ms->selectmoves(localchain, &moves);
			nloops = moves.nloops();
		}
		if (nloops == 0) {
			std::cout << "Cannot generate trial moves in BackBoneMinimizer."
					<< std::endl;
			return std::shared_ptr < Chain > (nullptr);
		}
		EnergyComponents pecomp0;
		temp = control.gettemperature(nstep_);
		for (int l = 0; l < nloops; ++l) {
			EnergyComponents pecomp1;
			double de = energy.deltaE(moves, l,&ce, &pecomp0, &pecomp1);
			double declash=pecomp1.energies[EnergyComponents::CLASH]-
					pecomp0.energies[EnergyComponents::CLASH];
			double paccept;
			if(declash > 1.0) paccept=0.0;
			else {
				paccept = exp(-(de-declash) / temp);
			}
			if (rng.realrng(0.0,1.0)()< paccept) {
				moves.updatechain(&localchain,ce.ssaspbtype(),l);
				moves.energy0 += de;
				ene += de;
				EnergyComponents &ecomp=ce.ecomp();
				for (int eidx = 0; eidx < ecomp.energies.size(); ++eidx) {
					ecomp.energies[eidx] += pecomp1.energies[eidx]
							- pecomp0.energies[eidx];
//					std::cout <<" "<<pecomp1.energies[eidx] <<":"<<pecomp0.energies[eidx];
				}
//				std::cout <<std::endl;
/*				pecomp0.energies = pecomp1.energies;
				EnergyComponents ecompsave=ce.ecomp();
				EnergyComponents ecomprecalc;
				ce.ecomp()=ecomprecalc;
				energy.totalenergy(localchain,&ce);
				int eidx=0;
				if(ecompsave.energies[0] - ce.energy(0) >1.0 ||
					ecompsave.energies[0]-ce.energy(0) <-1.0	)  {
					std::cout <<moves.startposi <<" "<< moves.endposi<<std::endl;
				}
				for (auto e : ecompsave.energies)
								std::cout << "\t" << e<<":"<<ce.energy(eidx++);
				std::cout << std::endl;
				ce.ecomp()=ecompsave;
*/
			}
			if (ene < emin_) {
				emin_ = ene;
				nochangesteps_ = 0;
				minconf_ = std::shared_ptr < Chain > (new Chain(localchain));
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

