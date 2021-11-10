/*
 * backbonelinker.cpp
 *
 *  Created on: 2016年11月23日
 *      Author: hyliu
 */

#include <backbone/backbonelinker.h>
#include <pdbstatistics/phipsidistr.h>
#include "dstl/randomengine.h"

using namespace NSPproteinrep;
using namespace NSPpdbstatistics;

BackBoneLinker::BackBoneLinker(unsigned int linkerlength,
		const std::vector<BackBoneSite> & fixedNend,
		const std::vector<BackBoneSite> & fixedCend, unsigned int maxgly,
		unsigned int maxpro) :
		NSPgeometry::FixEndsLoop<typename ChainTreeTopology::AtomKey>(nullptr), maxgly_(
				maxgly), maxpro_(maxpro) {
	linkerlength_ = linkerlength;
	headlength_ = fixedNend.size();
	taillength_ = fixedCend.size();
	length_ = linkerlength_ + headlength_ + taillength_;
	chaintopo_ = makeMainChainHeavyTopo(length_);
	topo_ = static_cast<Topology*>(chaintopo_.get());
	crd_ = std::shared_ptr < ChainTreeCrd
			> (new ChainTreeCrd(chaintopo_.get()));
	coord_ =
			static_cast<typename NSPgeometry::FixEndsLoop<AtomKey>::Coordinate *>(crd_.get());
	crd_->initWithIdealIntCrd();
	success_ = false;
//	std::vector<typename ChainTreeCrd::BackBoneTorsions> tor=linkertorsions();
//	for(auto & t:tor) std::cout <<t.toString() <<std::endl;
	rmsd_ = 100000;
	for (unsigned int i = 0; i < fixedNend.size(); ++i) {
		double deg = 3.14159265358979323846 / 180.0;
		double omiga = 180.0;
		if (i > 0)
			omiga = fixedNend[i - 1].omiga();
		crd_->setBackBoneTorsionAt(i,
				ChainTreeCrd::BackBoneTorsions(fixedNend[i].phi() * deg,
						fixedNend[i].psi() * deg, omiga * deg));
		AtomKey n_k = chaintopo_->genKey(i,
				BackBoneData::atomnames[BackBoneData::N]);
		AtomKey ca_k = chaintopo_->genKey(i,
				BackBoneData::atomnames[BackBoneData::CA]);
		AtomKey c_k = chaintopo_->genKey(i,
				BackBoneData::atomnames[BackBoneData::C]);
		AtomKey o_k = chaintopo_->genKey(i,
				BackBoneData::atomnames[BackBoneData::O]);
		NSPgeometry::XYZ ncrd = fixedNend[i].getcrd(BackBoneSite::NCRD);
		NSPgeometry::XYZ ccrd = fixedNend[i].getcrd(BackBoneSite::CCRD);
		NSPgeometry::XYZ cacrd = fixedNend[i].getcrd(BackBoneSite::CACRD);
		NSPgeometry::XYZ ocrd = fixedNend[i].getcrd(BackBoneSite::OCRD);
		addStartCrd(n_k, ncrd);
		addStartCrd(ca_k, cacrd);
		addStartCrd(c_k, ccrd);
		addStartCrd(o_k, ocrd);
	}
	copyStartCrdsTo(coord_);

	for (unsigned int j = 0; j < fixedCend.size(); ++j) {
		unsigned int i = j + linkerlength_ + headlength_;
		AtomKey n_k = chaintopo_->genKey(i,
				BackBoneData::atomnames[BackBoneData::N]);
		AtomKey ca_k = chaintopo_->genKey(i,
				BackBoneData::atomnames[BackBoneData::CA]);
		AtomKey c_k = chaintopo_->genKey(i,
				BackBoneData::atomnames[BackBoneData::C]);
		AtomKey o_k = chaintopo_->genKey(i,
				BackBoneData::atomnames[BackBoneData::O]);
		NSPgeometry::XYZ ncrd = fixedCend[j].getcrd(BackBoneSite::NCRD);
		NSPgeometry::XYZ ccrd = fixedCend[j].getcrd(BackBoneSite::CCRD);
		NSPgeometry::XYZ cacrd = fixedCend[j].getcrd(BackBoneSite::CACRD);
		NSPgeometry::XYZ ocrd = fixedCend[j].getcrd(BackBoneSite::OCRD);
		addEndCrd(n_k, ncrd);
		addEndCrd(ca_k, cacrd);
		addEndCrd(c_k, ccrd);
		addEndCrd(o_k, ocrd);
		double deg = 3.14159265358979323846 / 180.0;
		double omiga = 180.0;
		if (j > 0)
			omiga = fixedCend[j - 1].omiga();
		crd_->setBackBoneTorsionAt(i,
				ChainTreeCrd::BackBoneTorsions(fixedCend[j].phi() * deg,
						fixedCend[j].psi() * deg, omiga * deg));
	}
	for (unsigned int m = 0; m < linkerlength_; ++m) {
		unsigned int posi = m + headlength_;
		AtomKey ca_k = chaintopo_->genKey(posi,
				BackBoneData::atomnames[BackBoneData::CA]);
		AtomKey c_k = chaintopo_->genKey(posi,
				BackBoneData::atomnames[BackBoneData::C]);
		addRotatable(ca_k);
		addRotatable(c_k);
	}
	isgly_.resize(length_, false);
	ispro_.resize(length_, false);
}

void BackBoneLinker::generateTopLinkers(unsigned int ncopies,
		unsigned int ncandidates) {
	NSPdstl::RandomEngine<> &rneg = NSPdstl::RandomEngine<>::getinstance();
	rneg.setrealrng(0, 1);
	success_ = false;
	rmsd_ = 100000;
	int ntry = 0;
	std::priority_queue<ScoredLinker> queue;
	while (ntry < ncopies * ncandidates) {
		chooseglypro();
		bool initial = false;
		unsigned int ntryinit = 0;
		double tmin=10000.0;
		while (!initial) {
			if (ntryinit++ > 50000000) {
				std::cout << "Choosing random phi psi failed. Something wrong"
						<< std::endl;
				exit(1);
			}
			setRandomRotations(coord_, rneg.realrng());
			double t= phipsiscore();
			if(t <tmin) tmin=t;
			initial = t < 0.4*(double)(linkerlength())-1.4;
			if (ntryinit == (ntryinit/50000u)*50000u)
			std::cout << ntryinit << " init phipsiscore " <<tmin <<std::endl;
		}

		std::cout << "Ntryinit: " <<ntryinit << " init phipsiscore: " <<phipsiscore() <<std::endl;
		LinkerScorer phipsiscorer(this);
		success_ = MCclosure(coord_, 2000 * linkerlength_, 0.2, phipsiscorer,
				&rmsd_);
		double tscore = phipsiscore();
		std::cout << "phipsi score: " << tscore << std::endl;
		if (success_) {
			if (queue.size() >= ncopies) {
				if (tscore > queue.top().phipsienergy)
					continue;
			}
			ScoredLinker sl;
			sl.phipsienergy = tscore;
			sl.torsions = linkertorsions();
			sl.rmsd = rmsd_;
			queue.push(sl);
			if (queue.size() > ncopies)
				queue.pop();
		}
		++ntry;
	}
	if (queue.size() == ncopies)
		success_ = true;
	storedlinkers_.clear();
	while (!queue.empty()) {
		storedlinkers_.push_back(queue.top());
		queue.pop();
	}
}
bool BackBoneLinker::popTopLinker() {
	if (storedlinkers_.empty())
		return false;
	success_ = true;
	rmsd_ = storedlinkers_.back().rmsd;
	setlinkertorsions(storedlinkers_.back().torsions);
	std::cout << "Best linker generated with RMSD = " << rmsd_
			<< " and phipsi energy = " << storedlinkers_.back().phipsienergy
			<< std::endl;
	for (int posi = 0; posi < length_; ++posi) {
		std::cout << crd_->getBackBoneTorsions(posi).toString() << std::endl;
	}
	storedlinkers_.pop_back();
	return success_;
}
void BackBoneLinker::chooseglypro() {
	NSPdstl::RandomEngine<> &rneg = NSPdstl::RandomEngine<>::getinstance();
	const double PGLY { 0.024 };
	const double PPRO { 0.049 };
	rneg.setrealrng(0, 1);
	for (int i = headlength_; i < headlength_ + linkerlength_; ++i) {
		isgly_[i] = false;
		ispro_[i] = false;
		crd_->setOmega(i, 3.14159265);
		double rn = rneg.realrng()();
		if (rn < PGLY) {
			isgly_[i] = true;
		} else if (rn < PGLY + PPRO) {
			ispro_[i] = true;
			double q = rneg.realrng()();
			if (q > 0.90)
				crd_->setOmega(i, 0);
		}
	}
}

double BackBoneLinker::phipsiscore() {
	std::pair<double, double> ecoil { -2.42, 1.78 };
	std::pair<double, double> egly { -2.61, 1.97 };
	std::pair<double, double> eprepro { -2.27, 1.17 };
	double rad = 180.0 / 3.14159265;
	double ene = 0.0;
	double e_av = 0.0, e_fluc = 0.0;
	for (int posi = headlength_; posi < headlength_ + linkerlength_; ++posi) {
		double phi = crd_->getPhi(posi) * rad;
		double psi = crd_->getPsi(posi) * rad;

		const PhiPsiDistr *distr;
		if (isgly_[posi]) {
			distr = &(PhiPsiDistr::glydistr());
			e_av = ecoil.first;
			e_fluc = ecoil.second * ecoil.second;
		} else if (ispro_[posi + 1]) {
			distr = &(PhiPsiDistr::preprodistr());
			e_av = egly.first;
			e_fluc = egly.second * egly.second;

		} else {
			distr = &(PhiPsiDistr::coildistr());
			e_av = eprepro.first;
			e_fluc = eprepro.second * eprepro.second;
		}

		ene += distr->statisticalenergy(phi, psi)-e_av;
//		std::cout <<"(" << phi <<" "<<psi <<") " <<ene<<std::endl;
	}
	ene /= (double) (linkerlength_);
	return ene;
}
/*
 bool BackBoneLinker::getLinker(std::vector<BackBoneSite> *linker,
 unsigned int ncandidates,
 bool regenerate){
 if(regenerate || !success_)  generateLinker(ncandidates);
 if(!success_) return false;
 makealllinkersites(linker);
 return success_;
 }
 */
unsigned int BackBoneLinker::getLinkers(unsigned int ncopies,
		unsigned int ncandidates, typename std::vector<std::vector<BackBoneSite>>::iterator linkeriter) {
	generateTopLinkers(ncopies, ncandidates);
	unsigned int ngenerated = storedlinkers_.size();
	for (unsigned int n = 0; n < ngenerated; ++n) {
		if (popTopLinker())
			makealllinkersites(&(*(linkeriter++)));
		else
			return false;
	}
	return ngenerated;
}

void BackBoneLinker::makealllinkersites(
		std::vector<BackBoneSite> *linker) const {
	const std::map<AtomKey, NSPgeometry::XYZ> & crdmap = coord_->crdmap();
	for (unsigned int m = 0; m < linkerlength_; ++m) {
		unsigned int posi = m + headlength_;
		AtomKey n_k = chaintopo_->genKey(posi,
				BackBoneData::atomnames[BackBoneData::N]);
		AtomKey ca_k = chaintopo_->genKey(posi,
				BackBoneData::atomnames[BackBoneData::CA]);
		AtomKey c_k = chaintopo_->genKey(posi,
				BackBoneData::atomnames[BackBoneData::C]);
		AtomKey o_k = chaintopo_->genKey(posi,
				BackBoneData::atomnames[BackBoneData::O]);
		std::vector<NSPgeometry::XYZ> sitecrd;
		sitecrd.push_back(crdmap.at(n_k));
		sitecrd.push_back(crdmap.at(ca_k));
		sitecrd.push_back(crdmap.at(c_k));
		sitecrd.push_back(crdmap.at(o_k));
		linker->push_back(makelinkersite(posi, sitecrd));
	}
}

BackBoneSite BackBoneLinker::makelinkersite(unsigned int posi,
		const std::vector<NSPgeometry::XYZ> & sitecrd) const {
	BackBoneSite bs;
	bs.pdbid="000";
	bs.chainid = 'A';
	bs.resid = posi;
	bs.resseq = posi;
	bs.sscode = 'C';
	bs.resname = "GLY";
	double rad=180.0/3.14159265358979323846;
	double t=crd_->getPhi(posi)*rad;
	while(t >180.0) t-=360.0;
	while (t <-180.0) t+=360.0;
	bs.data_[BackBoneSite::PHI] = t;
	t=crd_->getPsi(posi)*rad;
	while(t >180.0) t-=360.0;
	while (t <-180.0) t+=360.0;
	bs.data_[BackBoneSite::PSI] = t ;
	if (posi < length_ - 1){
		t=crd_->getOmega(posi + 1)*rad;
		while(t >180.0) t-=360.0;
		while (t <-180.0) t+=360.0;
		bs.data_[BackBoneSite::OMIGA] = t;
	}
	bs.changecrd(sitecrd);
	return std::move(bs);
}

