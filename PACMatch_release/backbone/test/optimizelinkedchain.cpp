/*
 * optimizelinkedchain.cpp
 *
 *  Created on: 2017年5月4日
 *      Author: hyliu
 */

#include "backbone/mainchain.h"
#include "backbone/energyfunctions.h"
#include "backbone/closingloop.h"
#include "dstl/randomengine.h"
#include "dstl/topn.h"
#include <algorithm>
#include <memory>
using namespace NSPproteinrep;
class AnnealingTemperature {
public:
	AnnealingTemperature(double t_high, double t_low, int stepmax) :
			t_high_(t_high), t_low_(t_low), nstepmax_(stepmax) {
		;
	}
	double temperature_cos(int step) {
		double stepi = step % nstepmax_;
		stepi=sqrt(stepi / (double) nstepmax_);
		double phi = stepi  * 3.14159265;
		return t_low_ + 0.5 * (t_high_ - t_low_) * (1 + cos(phi));
	}
	double temperature_staged(int step) {
		double stepi = (double) (step % nstepmax_)/(double) nstepmax_;
		for(int i=0; i<times_.size();++i){
			if(stepi < times_[i]) return t_low_+temps_[i]*(t_high_-t_low_);
		}
	}
private:
	double t_high_;
	double t_low_;
	int nstepmax_;
	std::vector<double> times_{0.1,0.2,0.9,1.01};
	std::vector<double> temps_{1.0,0.2,0.1,0.0};
};
void savetopn(NSPdstl::TopN<std::shared_ptr<MainChain>>&topn) {
	std::vector<std::pair<std::shared_ptr<MainChain>, double>> topres =
			NSPdstl::topN2vector(topn);
	int resnum = 0;
	for (auto & res : topres) {
		std::map<std::string, double> etmp;
		std::cout << "Lowest energy " << resnum << ": " << res.second <<std::endl;
		//				NSPproteinrep::calctotalenergy(*(res.first), erank, etmp,
		//							erank.emaxallowed) << std::endl;
		std::string outname = "optichain" + std::to_string(resnum) + ".pdb";
		std::ofstream ofs;
		ofs.open(outname.c_str());
		writeSitesToPDB(ofs, *(res.first));
		ofs.close();
		std:outname = "optichain" + std::to_string(resnum++) + ".dat";
		ofs.open(outname.c_str());
		res.first->write(ofs);
	}
}
struct LowestEnergy {
	double step;
	double energy;
	bool isnewlowest(double e) {
		++step;
		if (e < energy) {
			energy = e;
			return true;
		} else
			return false;
	}
};
int main(int argc, char **argv) {
	EnergyTerms erank;
	EnergyTerms escreen;
	escreen.addstericclash();
	escreen.addtorsionene();
//	erank.addstericclash(100.0);
//	erank.addtorsionvecene();
	erank.addtetrasefene(0.17);
	NSPdstl::RandomEngine<> & rneg = NSPdstl::RandomEngine<>::getinstance();
	MainChain chain;
	chain.read(std::string(argv[1]));
	ChangeRegionSelector crs(&chain);
	double t_high = std::stod(std::string(argv[2]));
	double t_low = std::stod(std::string(argv[3]));
	int annealingsteps = std::stoi(std::string(argv[4]));
	int maxstep = std::stoi(std::string(argv[5]));
	AnnealingTemperature ant(t_high, t_low, annealingsteps);

//	double eclash = NSPproteinrep::calctotalenergy(chain, escreen,
//			total_energies, escreen.emaxallowed);
//	if (eclash > escreen.emaxallowed) {
//		std::cout << "Initial structure is not clash free." << std::endl;
//		exit(1);
//	}
	int loopstart = crs.flexiblebegin();
	int looplength = crs.flexibleend() - loopstart + 1;
	std::map<std::string, double> total_energies;
	double esave = calcloopenergy(chain, loopstart, looplength,
			std::vector<BackBoneSite>(), erank, total_energies,
			erank.emaxallowed);
	for(auto &e:total_energies){
		std::cout <<e.first<< "\t" <<e.second<<std::endl;
	}
//	double esave = NSPproteinrep::calctotalenergy(chain, erank, total_energies,
//			erank.emaxallowed);
	int nsave = 50;
	NSPdstl::TopN<std::shared_ptr<MainChain>> topn(nsave);
//	int maxstep = 400;
	int step = 0;
	std::cout << step << "\t" << esave << std::endl;
	LowestEnergy le { 0, 1000000000.0 };
	std::shared_ptr<MainChain> tosave;
	while (step < maxstep) {
		++step;
		if ((step % 500) == 0)
			savetopn(topn);
		int changestart, changeend;
		if (!crs.selectinoneloop(&changestart, &changeend))
			continue;
		double temper = ant.temperature_cos(step);
		int clmode=ClosingLoop::KEEPSEQ_MUTATECONF;;
		if(temper>1.0) clmode= ClosingLoop::KEEPSEQ_NEWCONF;
		ClosingLoop cl(chain, changestart, changeend,
				clmode);
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
		if (screenedsol.empty())
			continue;

		std::map<std::string, double> old_loop_energies;
		std::vector<BackBoneSite> tmp(changeend - changestart, BackBoneSite());
		std::copy(chain.begin() + changestart, chain.begin() + changeend,
				tmp.begin());
		double eold = calcloopenergy(chain, changestart,
				changeend - changestart, tmp, erank, old_loop_energies,
				erank.emaxallowed);
//		std::cout <<"p1: "<< esave <<"\t" <<eold<<std::endl;
		esave = esave - eold;
		for(auto & e:total_energies) {
			e.second = e.second-old_loop_energies[e.first];
		}
		int acceptedis = -1;
		std::shared_ptr<std::vector<BackBoneSite>> acceptedloop;
		for (int is = 0; is < screenedsol.size(); ++is) {
			std::shared_ptr<std::vector<BackBoneSite>> newloopi =
					screenedsol[is];
			std::map<std::string, double> new_loop_energies;
			double enew = calcloopenergy(chain, changestart,
					changeend - changestart, *newloopi, erank,
					new_loop_energies, erank.emaxallowed);
			if (enew > erank.emaxallowed)
				continue;
			double de = exp((eold - enew) / temper);
			double r = rneg.realrng(0, 1.0)();
			if (r < de) {
				acceptedis = is;
				acceptedloop = newloopi;
				eold = enew;
				old_loop_energies=new_loop_energies;
			}
		}
		esave += eold;
		for(auto & e:total_energies) {
			e.second = e.second+old_loop_energies[e.first];
		}
//		std::cout <<"p2: "<< esave <<"\t" <<eold<<std::endl;
		if (acceptedis >= 0) {
			std::copy(acceptedloop->begin(), acceptedloop->end(),
					chain.begin() + changestart);
			chain.resetresseq();
//			std::map<std::string,double> etmp;
//			double etot=NSPproteinrep::calctotalenergy(chain, erank, etmp,
//									erank.emaxallowed);
//			std::cout <<"Etot: "<< etot <<std::endl;
			if (le.isnewlowest(esave)) {
				tosave = std::shared_ptr < MainChain > (new MainChain(chain));
			}
			if (le.step >= 50) {
				if (topn.keep(le.energy)) {
					topn.push(tosave, le.energy);
					le.step = 0;
					if (topn.size() < nsave)
						le.energy = 100000000.0;
					else
						topn.top(&le.energy);
				}
			}
		}
		std::cout << step << "\t" << esave	<< std::endl;
		for(auto &e:total_energies){
			std::cout <<e.first<< "\t" <<e.second<<std::endl;
		}
	}
	savetopn(topn);
}

