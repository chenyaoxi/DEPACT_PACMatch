/*
 * chainloop.h
 *
 *  Created on: 2016年11月23日
 *      Author: hyliu
 */

#ifndef BACKBONE_BACKBONELINKER_H_
#define BACKBONE_BACKBONELINKER_H_

#include <backbone/backbonesite.h>
#include "proteinrep/chaintree.h"
#include "geometry/fixendsloop.h"
#include <queue>
namespace NSPproteinrep {

class BackBoneLinker: public NSPgeometry::FixEndsLoop<
		typename ChainTreeTopology::AtomKey> {
public:
	typedef typename ChainTreeTopology::AtomKey AtomKey;
	struct ScoredLinker {
		std::vector<ChainTreeCrd::BackBoneTorsions> torsions;
		double phipsienergy;
		double rmsd;
		bool operator<(const ScoredLinker & sl2) const {
			return phipsienergy <sl2.phipsienergy;
		}
	};
/*
	struct ScoredLinkerOrder {
		bool operator()(const ScoredLinker &sl1, const ScoredLinker &sl2) {
			return sl1.phipsienergy < sl2.phipsienergy;
		}
	};
	*/
	BackBoneLinker() :
			NSPgeometry::FixEndsLoop<typename ChainTreeTopology::AtomKey>(
					nullptr) {
		;
	}
	BackBoneLinker(unsigned int linkerlength,
			const std::vector<BackBoneSite> & fixedNend,
			const std::vector<BackBoneSite> &fixedCend, unsigned int maxgly = 0,
			unsigned int maxpro = 0);
//	void generateLinker(unsigned int ncandidates);
	void generateTopLinkers(unsigned int ncopies,
			unsigned int ncandidates);
	bool popTopLinker();
//	bool getLinker(std::vector<BackBoneSite> *linker,unsigned int ncandidates,
//			bool regenerate=false);
	unsigned int getLinkers(unsigned int ncopies, unsigned int ncandidates,
			typename std::vector<std::vector<BackBoneSite>>::iterator linkeriter);
	bool success() const {
		return success_;
	}
	double rmsd() const {
		return rmsd_;
	}
	void chooseglypro();
	double phipsiscore();
	ChainTreeCrd *crd() const {
		return crd_.get();
	}
	unsigned int linkerlength() const {
		return linkerlength_;
	}
	std::vector<ChainTreeCrd::BackBoneTorsions> linkertorsions() const {
		std::vector<ChainTreeCrd::BackBoneTorsions> torsions;
		for (unsigned int posi = headlength_;
				posi < headlength_ + linkerlength_; ++posi) {
			torsions.push_back(crd_->getBackBoneTorsions(posi));
		}
		return std::move(torsions);
	}
	void setlinkertorsions(
			const std::vector<ChainTreeCrd::BackBoneTorsions> &torsions) {
		for (unsigned int posi = 0; posi < linkerlength_; ++posi) {
			crd_->setBackBoneTorsionAt(posi + headlength_, torsions[posi]);
			crd_->calcXYZ();
		}
	}
private:
	std::shared_ptr<ChainTreeTopology> chaintopo_;
	unsigned int headlength_;
	unsigned int linkerlength_;
	unsigned int taillength_;
	unsigned int length_;
	unsigned int maxgly_;
	unsigned int maxpro_;
	std::shared_ptr<ChainTreeCrd> crd_;
	NSPgeometry::FixEndsLoop<AtomKey>::Coordinate *coord_;
	bool success_ { false };
	double rmsd_ { 1000000.0 };
	unsigned int buffersize_ { 200 };
	std::vector<ScoredLinker> storedlinkers_;
	std::vector<bool> isgly_;
	std::vector<bool> ispro_;
	BackBoneSite makelinkersite(unsigned int posi,
			const std::vector<NSPgeometry::XYZ> &crd) const;
	void makealllinkersites(std::vector<BackBoneSite> *linker) const;
};

class LinkerScorer {
public:
	LinkerScorer(BackBoneLinker *bl) :
			bl_(bl) {
		;
	}
	double operator()(const typename BackBoneLinker::AtomKey & ka, double rot) {
		double scale=4.0;
		bl_->crd()->rotateBond(ka, rot);
		double score = (bl_->phipsiscore()) * (double) (bl_->linkerlength());
		bl_->crd()->rotateBond(ka, -rot);
		return scale*score;
	}
private:
	BackBoneLinker *bl_;
};

}

#endif /* BACKBONE_BACKBONELINKER_H_ */
