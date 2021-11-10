/*
 * paretoloopopt.cpp
 *
 *  Created on: 2017年5月8日
 *      Author: hyliu
 */

#include "backbone/mainchain.h"
#include "backbone/energyfunctions.h"
#include "backbone/closingloop.h"
#include "dstl/randomengine.h"
#include "dstl/topn.h"
#include "dstl/paretofront.h"
#include <algorithm>
#include <memory>
using namespace NSPproteinrep;
void printfront(NSPdstl::FrontSets<std::shared_ptr<MainChain>> &front) {
	std::set<NSPdstl::FrontSets<std::shared_ptr<MainChain>>::MultiScoredType> s =
			front.saved();
	std::cout << "Energies of saved chains: " << std::endl;
	int id = 0;
	for (auto &c : s) {
		std::cout << id;
		double etot = 0.0;
		for (auto ene : c.scores_) {
			etot += ene;
			std::cout << "\t" << ene;
		}
		std::cout << "\t" << etot << std::endl;
		std::string output = "front" + std::to_string(id) + ".pdb";
		std::ofstream ofs;
		ofs.open(output.c_str());
		writeSitesToPDB(ofs, *(c.obj_));
		ofs.close();
		output = "front" + std::to_string(id) + ".dat";
		ofs.open(output.c_str());
		c.obj_->write(ofs);
		ofs.close();
		++id;
	}
}
void picknewstart(NSPdstl::FrontSets<std::shared_ptr<MainChain>> &front,
		MainChain *chain, std::map<std::string, double> *enes) {
	int i =
			NSPdstl::RandomEngine<>::getinstance().intrng(0, front.nsaved() - 1)();
	const NSPdstl::FrontSets<std::shared_ptr<MainChain>>::MultiScoredType & newstart =
			front.getsaved(i);
	*chain = *(newstart.obj_);
	enes->at("TetraSefEne") = newstart.scores_[0];
	enes->at("TorsionVecEne") = newstart.scores_[1];
}

int main(int argc, char **argv) {
	EnergyTerms erank;
	EnergyTerms escreen;
	escreen.addstericclash(erank.emaxallowed + 100.0);
	escreen.addtorsionene();
//	erank.addstericclash(100.0);
	erank.addtorsionvecene();
	erank.addtetrasefene(0.17);
	NSPdstl::RandomEngine<> & rneg = NSPdstl::RandomEngine<>::getinstance();
	MainChain chain;
	chain.read(std::string(argv[1]));
	ChangeRegionSelector crs(&chain);
	int ene_dim = 2;
	NSPdstl::FrontSets<std::shared_ptr<MainChain>> front(ene_dim);
	int ndir = std::stod(std::string(argv[2]));
	int ntop = std::stod(std::string(argv[3]));
	front.generatedirections2D(ndir, ntop);
	int maxstep = std::stoi(std::string(argv[4]));
	int seed = std::stoi(std::string(argv[5]));
	rneg.init(seed);
	int loopstart = crs.flexiblebegin();
	int looplength = crs.flexibleend() - loopstart + 1;
	std::map<std::string, double> total_energies;
	double esave = NSPproteinrep::calctotalenergy(chain, erank, total_energies,
			erank.emaxallowed);
//	double esave = calcloopenergy(chain, loopstart, looplength,
//			std::vector<BackBoneSite>(), erank, total_energies,
//			erank.emaxallowed);
	for (auto &e : total_energies) {
		std::cout << e.first << "\t" << e.second << std::endl;
	}

	int nsave = 10;
	NSPdstl::TopN<std::shared_ptr<MainChain>> topn(nsave);
//	int maxstep = 400;
	int step = 0;
	std::cout << step << "\t" << esave << std::endl;
	std::shared_ptr<MainChain> tosave;
	int nblocked = 0;
	double runningav_nblock=0.0;
	while (step < maxstep) {
		++step;
		if (step % 100 == 0)
			printfront(front);
		int changestart, changeend;
		if (!crs.selectinoneloop(&changestart, &changeend))
			continue;
		int clmode = ClosingLoop::MUTATESEQ_CONF;
		;
		double r = NSPdstl::RandomEngine<>::getinstance().realrng(0, 1.0)();
		if (r < 0.5)
			clmode = ClosingLoop::NEWSEQ_NEWCONF;
		ClosingLoop cl(chain, changestart, changeend, clmode);
		int nscreenedsol = 0;
		int ntry = 0;
		std::vector<std::shared_ptr<std::vector<BackBoneSite>>>screenedsol;

		while (screenedsol.empty() && ntry < 100) {
			++ntry;
			int nsol = cl.solve();
			for (int i = 0; i < nsol; ++i) {
				std::shared_ptr<std::vector<BackBoneSite>> newloopi(
						new std::vector<BackBoneSite>());
				cl.getsitessolution(i, newloopi.get());
				std::map<std::string, double> energies;
				double escr = calcloopenergy(chain, changestart,
						newloopi->size(), *newloopi, escreen, energies,
						escreen.emaxallowed);
				if (escr >= escreen.emaxallowed
						|| energies["TorsionEne"]
								> -1.0 * (double) (newloopi->size()))
					continue;
				screenedsol.push_back(newloopi);
			}
		}
		if (screenedsol.empty()) {
			continue;
		}
		std::cout << step << "\t" << screenedsol.size() << std::endl;
		std::map<std::string, double> old_loop_energies;
		std::vector<BackBoneSite> tmp(changeend - changestart, BackBoneSite());
		std::copy(chain.begin() + changestart, chain.begin() + changeend,
				tmp.begin());
		double eold_p = calcloopenergy(chain, changestart,
				changeend - changestart, tmp, erank, old_loop_energies,
				erank.emaxallowed);
//		std::cout <<"p1: "<< esave <<"\t" <<eold<<std::endl;
		double eold = eold_p;
		esave = esave - eold;
		for (auto & e : total_energies) {
			e.second = e.second - old_loop_energies[e.first];
		}
		bool moved = false;
		for (int is = 0; is < screenedsol.size(); ++is) {
			std::shared_ptr<std::vector<BackBoneSite>> newloopi =
					screenedsol[is];
			std::map<std::string, double> new_loop_energies;
			double enew = calcloopenergy(chain, changestart,
					changeend - changestart, *newloopi, erank,
					new_loop_energies, erank.emaxallowed);
			if (enew > erank.emaxallowed)
				continue;
			if (enew - eold_p > -1.e-5 && enew - eold_p < 1.e-5)
				continue;
			std::vector<double> ene_new;
			ene_new.push_back(
					total_energies["TetraSefEne"]
							+ new_loop_energies["TetraSefEne"]);
			ene_new.push_back(
					total_energies["TorsionVecEne"]
							+ new_loop_energies["TorsionVecEne"]);
//			NSPdstl::FrontSets<std::shared_ptr<MainChain>>::FrontSet *fs;
//			double cs;
//			if (front.willsavetodir(ene_new, fs, &cs)) {
			if (front.willsave(ene_new)) {
				moved = true;
				std::copy(newloopi->begin(), newloopi->end(),
						chain.begin() + changestart);
				chain.resetresseq();
				tosave = std::shared_ptr < MainChain > (new MainChain(chain));
//				front.savetodir(tosave, ene_new, *fs, cs);
				front.save(tosave,ene_new);
				eold = enew;
				old_loop_energies = new_loop_energies;
				std::cout << step << "\t" << esave + eold << std::endl;
				for (auto &e : total_energies) {
					std::cout << e.first << "\t"
							<< e.second + old_loop_energies[e.first]
							<< std::endl;
				}
			}
		}
		esave = esave + eold;
		for (auto & e : total_energies) {
			e.second = e.second + old_loop_energies[e.first];
		}
		if (moved) {
			runningav_nblock= 0.8*runningav_nblock + 0.2*(double)nblocked;
			nblocked = 0;
		} else {
			++nblocked;
			if (nblocked >= 20 ) {
				runningav_nblock= 0.8*runningav_nblock + 0.2*(double)nblocked;
				std::cout << "Start from another configuration in stored front. "
						<< "avnblock: " << runningav_nblock <<std::endl;
				picknewstart(front, &chain, &total_energies);
				esave = total_energies["TetraSefEne"]
						+ total_energies["TorsionVecEne"];
							nblocked=0;
			}
		}
		if(runningav_nblock > 16) {
			printfront(front);
			break;
		}
	}
}

