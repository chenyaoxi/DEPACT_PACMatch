/*
 * closealoop.h
 *
 *  Created on: 2017年8月1日
 *      Author: hyliu
 */

#ifndef BACKBONE_CLOSEALOOP_H_
#define BACKBONE_CLOSEALOOP_H_
#include "backbone/backbonesite.h"
#include "loopclosure/loopsolver.h"
#include "pdbstatistics/phipsidistr.h"
#include <memory>
#include <set>
namespace NSPproteinrep {
typedef std::vector<BackBoneSite> Loop;
class CloseALoop {
public:
	CloseALoop() {
		;
	}
	CloseALoop(const std::vector<BackBoneSite> & hostchain, int loopstart,
			int looplength, const std::vector<double> &inittorsions);
	CloseALoop(const std::vector<BackBoneSite> & hostchain, int loopstart,
			int looplength);
	std::vector<std::shared_ptr<Loop>> getsolutions(
			const std::set<int> &fixedpositions = std::set<int>());
private:
	const std::vector<BackBoneSite> *hostchain_ { nullptr };
	int loopstart_ { -1 };
	int looplength_ { -1 };
	Loop refloop_;
	std::vector<NSPloopclosure::LoopSolution> solutions_;
	std::vector<int> p2ofsolutions_;
	std::vector<NSPgeometry::XYZ> fixcrds_;
	int findsolutions(int posi);
	void buildasolution(int idx, Loop *result);
	void init(const std::vector<double> &inittorsions);
};

/*!Linking two given "ending" backbone sites with peptide backbones
 * comprising a series of "building blocks".
 * A building block can be a user-supplied one with fixed internal phi,psi and omiga angles.
 * it can be also a loop with random phi,psi and trans-omiga angles sampled
 * in methods of this class.
 * Each rigid building block must be flanked by random loops at both sides, while each random loop must be flanked
 * by rigid blocks or ending sites at both sides.
 *
 * when solving the loop closure problem, the phi,psi and omiga angles of the ending sites and
 * sites of the rigid building blocks are not changed. The coordinates of the
 * N and C_alpha atoms of the first linker site are computed using backbone atom coordinates of the left or N-terminal
 * ending site together with it psi and
 * omiga angles. The coordinates of the Calpha and C atoms of the last linker site are computed
 * using the backbone atom coordinates of the right or C-terminal ending site together with its phi
 * angles.
 * Within the random blocks, the peptide bonds are generated in the trans configuration.
 * The phi,psi angles are sampled from the statistical PDB distributions, including
 * distributions coming from GLY and PRO.
 *
 */
class LinkWithBlocks {
public:
	struct BackBoneTorsions {
		double phi;
		double psi;
		double omiga;
		BackBoneTorsions() {
			;
		}
		BackBoneTorsions(double p, double ps, double o = 180.0) :
				phi(p), psi(ps), omiga(o) {
			;
		}
	};
	LinkWithBlocks(const BackBoneSite &leftsite, const BackBoneSite &rightsite) :
			left_(leftsite), right_(rightsite) {
		;
	}
	void reset() {
		blocklengths_.clear();
		fixedtorsions_.clear();
		linkers_.clear();
	}
	;
	void addcoilblock(int blocklength) {
		blocklengths_.push_back(blocklength);
	}
	void addrigidblock(const std::vector<BackBoneTorsions> &torsions) {
		int len = torsions.size();
		int id = blocklengths_.size();
		blocklengths_.push_back(len);
		fixedtorsions_.insert(std::make_pair(id, torsions));
	}
	void addrigidblock(const std::vector<BackBoneSite> segment) {
		std::vector<BackBoneTorsions> torsions;
		for (auto &s : segment) {
			torsions.push_back(BackBoneTorsions(s.phi(), s.psi(), s.omiga()));
		}
		addrigidblock(torsions);
	}
	template<typename RNG>
	void addhelixblock(RNG &rng, int length){
		std::vector<BackBoneTorsions> torsions;
		NSPpdbstatistics::PhiPsiDistr & distr =
				NSPpdbstatistics::PhiPsiDistr::helixdistr();
		for(int i=0;i<length;++i){
			double phi,psi;
			distr.randomphipsi(rng,&phi,&psi);
			torsions.push_back(BackBoneTorsions(phi,psi));
		}
		addrigidblock(torsions);
	}
	template<typename RNG>
	int findnewlinkers(RNG &rng) {
		int looplength = 0;
		for (auto l : blocklengths_)
			looplength += l;
		assert(looplength >= 3);
		std::vector<BackBoneSite> hostchain(looplength + 2);
		hostchain[0] = left_;
		hostchain.back() = right_;
		std::vector<double> inittorsions = sample_torsions(rng);
		int tidx = inittorsions.size() - 3;
		genprevbackbonesite(&(hostchain[looplength + 1]),
				inittorsions[tidx + 2], inittorsions[tidx + 1],
				inittorsions[tidx], &(hostchain[looplength]));
		CloseALoop closer(hostchain, 1, looplength, inittorsions);
		std::vector<std::shared_ptr<Loop>> sols = closer.getsolutions();
		for (auto &s : sols) {
			linkers_.push_back(s);
		}
		return sols.size();
	}

	int nlinkers() const {
		return linkers_.size();
	}

	template<typename RNG>
	std::shared_ptr<Loop> poplinker(RNG &rng, int ntry = 1) {
		if (linkers_.empty()) {
			while (ntry > 0) {
				--ntry;
				int nl = findnewlinkers(rng);
				if (nl > 0)
					break;
			}
			if (linkers_.empty())
				return nullptr;
		}
		std::shared_ptr<Loop> res = linkers_.back();
		linkers_.pop_back();
		return res;
	}
private:
	BackBoneSite left_;
	BackBoneSite right_;
	std::vector<int> blocklengths_; //the lengths of blocks forming the linker,including
									//both random coils and user-specified rigid blocks
	std::map<int, std::vector<BackBoneTorsions>> fixedtorsions_; //The keys are the sequential block ids
	//of the rigid blocks. The values are
	//sequential lists of backbone torsions
	// for individual rigid blocks.
	std::vector<std::shared_ptr<Loop>> linkers_;
	template<typename RNG>
	std::vector<double> sample_torsions(RNG &rng) {
		std::vector<double> torsions;
		NSPpdbstatistics::PhiPsiDistr & distr =
				NSPpdbstatistics::PhiPsiDistr::mixcoildistr();
		double omiga_trans = 180.0;
		for (int i = 0; i < blocklengths_.size(); ++i) {
			if (fixedtorsions_.find(i) != fixedtorsions_.end()) {
				auto & ft = fixedtorsions_.at(i);
				for (auto &bt : ft) {
					torsions.push_back(bt.phi);
					torsions.push_back(bt.psi);
					torsions.push_back(bt.omiga);
				}
			} else {
				for (int j = 0; j < blocklengths_[i]; ++j) {
					double phi, psi;
					distr.randomphipsi(rng, &phi, &psi);
					torsions.push_back(phi);
					torsions.push_back(psi);
					torsions.push_back(omiga_trans);
				}
			}
		}
		return torsions;
	}
	;
};
}

#endif /* BACKBONE_CLOSEALOOP_H_ */
