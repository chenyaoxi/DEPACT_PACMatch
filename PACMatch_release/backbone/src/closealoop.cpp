/*
 * closealoop.cpp
 *
 *  Created on: 2017年8月1日
 *      Author: hyliu
 */

#include "backbone/closealoop.h"
#include "pdbstatistics/phipsidistr.h"
#include "dstl/randomengine.h"
#include <algorithm>
using namespace NSPproteinrep;
using namespace NSPgeometry;

CloseALoop::CloseALoop(const std::vector<BackBoneSite> & hostchain,
		int loopstart, int looplength, const std::vector<double> &inittorsions) :
		hostchain_(&hostchain), loopstart_ { loopstart }, looplength_(
				looplength) {
	init(inittorsions);
}
CloseALoop::CloseALoop(const std::vector<BackBoneSite> & hostchain,
		int loopstart, int looplength) :
		hostchain_(&hostchain), loopstart_ { loopstart }, looplength_(
				looplength) {
	std::vector<double> inittorsions;
	NSPdstl::RandomEngine<> &rneg = NSPdstl::RandomEngine<>::getinstance();
	rneg.setrealrng(0,1);
	for (int i = 0; i < looplength_; ++i) {
		const NSPpdbstatistics::PhiPsiDistr *distr =
				&(NSPpdbstatistics::PhiPsiDistr::coildistr());
		std::string resname = hostchain[loopstart + i].resname;
		if (resname == "GLY")
			distr = &(NSPpdbstatistics::PhiPsiDistr::glydistr());
		else if (resname == "PRO")
			distr = &(NSPpdbstatistics::PhiPsiDistr::transprodistr());
		else if (hostchain[loopstart + i + 1].resname == "PRO")
			distr = &(NSPpdbstatistics::PhiPsiDistr::preprodistr());
		double phi, psi;
		distr->randomphipsi(rneg.realrng(), &phi, &psi);
		if(phi>180.0) phi -=360.0;
		if(psi>180.0) psi -=360.0;
		inittorsions.push_back(phi);
		inittorsions.push_back(psi);
		inittorsions.push_back(180.0);
	}
	init(inittorsions);
}

void CloseALoop::init(const std::vector<double> &inittorsions) {
	assert(loopstart_ > 0);
	assert(looplength_ >= 3);
	assert(loopstart_ + looplength_ < hostchain_->size());
	refloop_.resize(looplength_, BackBoneSite());
	BackBoneSite site0((*hostchain_)[loopstart_ - 1]);
	BackBoneSite *psite = &site0;
	bool cispep = site0.nextpepcis();
	int tidx = 0;
	for (int i = 0; i < looplength_; ++i) {
		BackBoneSite *nsite = &(refloop_[i]);
		genbackbonesite(psite, cispep, inittorsions[tidx],
				inittorsions[tidx + 1], nsite);
		nsite->resname = (*hostchain_)[loopstart_ + i].resname;
		psite = nsite;
		tidx += 2;
		double omiga = inittorsions[tidx];
		cispep = omiga > -90.0 && omiga < 90.0;
		++tidx;
	}
	refloop_.back().data_[BackBoneSite::OMIGA] = (*hostchain_)[loopstart_
			+ looplength_].omiga();
	fixcrds_.clear();
	fixcrds_.push_back((*hostchain_)[loopstart_].ncrd());
	fixcrds_.push_back((*hostchain_)[loopstart_].cacrd());
	fixcrds_.push_back((*hostchain_)[loopstart_ + looplength_ - 1].cacrd());
	fixcrds_.push_back((*hostchain_)[loopstart_ + looplength_ - 1].ccrd());
	solutions_.clear();
	p2ofsolutions_.clear();
}
std::vector<std::shared_ptr<Loop>> CloseALoop::getsolutions(
		const std::set<int> & fixedpositions) {
//	solutions_.clear();
//	p2ofsolutions_.clear();
	int nsol_old = 0;
	for (int posi = 1; posi < looplength_ - 1; posi++) {
		if (fixedpositions.find(posi + loopstart_) != fixedpositions.end())
			continue;
		int nsol = findsolutions(posi);
		for (int n = nsol_old; n < nsol; ++n)
			p2ofsolutions_.push_back(posi);
		nsol_old = nsol;
	}
	std::vector<std::shared_ptr<Loop>> results;
	for (int n = 0; n < solutions_.size(); ++n) {
		results.push_back(std::shared_ptr < Loop > (new Loop));
		buildasolution(n, results.back().get());
	}
	return results;
}
void CloseALoop::buildasolution(int idx, Loop *result) {
	NSPloopclosure::LoopSolution &sol = solutions_.at(idx);
	int p2 = p2ofsolutions_.at(idx);
	const NSPgeometry::RigidTransform & rt1 = sol.rt1;
	const NSPgeometry::RigidTransform & rt2 = sol.rt2;
	int length = refloop_.size();
	result->resize(looplength_, BackBoneSite());
	std::copy(refloop_.begin(), refloop_.end(), result->begin());
	for (int p = 0; p < length; ++p) {
		std::vector<NSPgeometry::XYZ> crd;
		if (p == 0) {
			crd.push_back(fixcrds_[0]);
			crd.push_back(fixcrds_[1]);
			crd.push_back(rt1.applytoCopy(refloop_[p].ccrd()));
			crd.push_back(rt1.applytoCopy(refloop_[p].ocrd()));
		} else if (p < p2) {
			crd.push_back(rt1.applytoCopy(refloop_[p].ncrd()));
			crd.push_back(rt1.applytoCopy(refloop_[p].cacrd()));
			crd.push_back(rt1.applytoCopy(refloop_[p].ccrd()));
			crd.push_back(rt1.applytoCopy(refloop_[p].ocrd()));
		} else if (p == p2) {
			double rca_c1 = NSPgeometry::distance(refloop_[p].cacrd(),
					refloop_[p].ccrd());
			crd.push_back(rt1.applytoCopy(refloop_[p].ncrd()));
			crd.push_back(rt1.applytoCopy(refloop_[p].cacrd()));
//
			crd.push_back(rt2.applytoCopy(refloop_[p].ccrd()));
			crd.push_back(rt2.applytoCopy(refloop_[p].ocrd()));
			double rca_c2 = NSPgeometry::distance(
					rt1.applytoCopy(refloop_[p].cacrd()),
					rt2.applytoCopy(refloop_[p].ccrd()));
			if (rca_c2 > 1.6) {
				std::cout << "rca_c at p2: " << rca_c1 << "\t" << rca_c2
						<< std::endl;
			}
		} else if (p != length - 1) {
			crd.push_back(rt2.applytoCopy(refloop_[p].ncrd()));
			crd.push_back(rt2.applytoCopy(refloop_[p].cacrd()));
			crd.push_back(rt2.applytoCopy(refloop_[p].ccrd()));
			crd.push_back(rt2.applytoCopy(refloop_[p].ocrd()));
		} else {
			crd.push_back(rt2.applytoCopy(refloop_[p].ncrd()));
			crd.push_back(fixcrds_[2]);
			crd.push_back(fixcrds_[3]);
			XYZ rn4 = hostchain_->at(loopstart_ + looplength_).ncrd();
			NSPgeometry::XYZ ocrd = NSPgeometry::InternaltoXYZ(fixcrds_[3],
					fixcrds_[2], rn4, 1.23, 2.164, 3.14159265);
			crd.push_back(ocrd);
		}
		result->at(p).changecrd(crd.begin());
		result->at(p).sscode = 'C';
	}
	double phi1, psi1, phi2, psi2, phi3, psi3;
	double rad = 180.0 / 3.14159265;
	BackBoneSite &s1 = result->at(0);
	XYZ rc0 = hostchain_->at(loopstart_ - 1).ccrd();
	phi1 = NSPgeometry::torsion(rc0, s1.ncrd(), s1.cacrd(), s1.ccrd());
	psi1 = NSPgeometry::torsion(s1.ncrd(), s1.cacrd(), s1.ccrd(),
			result->at(1).ncrd());
	s1.data_[BackBoneSite::PHI] = phi1 * rad;
	s1.data_[BackBoneSite::PSI] = psi1 * rad;
	BackBoneSite &s2 = result->at(p2);
	phi2 = NSPgeometry::torsion(result->at(p2 - 1).ccrd(), s2.ncrd(),
			s2.cacrd(), s2.ccrd());
	psi2 = NSPgeometry::torsion(s2.ncrd(), s2.cacrd(), s2.ccrd(),
			result->at(p2 + 1).ncrd());
	s2.data_[BackBoneSite::PHI] = phi2 * rad;
	s2.data_[BackBoneSite::PSI] = psi2 * rad;
	BackBoneSite &s3 = result->back();
	phi3 = NSPgeometry::torsion(result->at(length - 2).ccrd(), s3.ncrd(),
			s3.cacrd(), s3.ccrd());
	XYZ rn4 = hostchain_->at(loopstart_ + looplength_).ncrd();
	psi3 = NSPgeometry::torsion(s3.ncrd(), s3.cacrd(), s3.ccrd(), rn4);
	s3.data_[BackBoneSite::PHI] = phi3 * rad;
	s3.data_[BackBoneSite::PSI] = psi3 * rad;
}
int CloseALoop::findsolutions(int posi) {
	assert(posi > 0 && posi < looplength_ - 1);
	std::vector<NSPgeometry::XYZ> atomcrds;
	atomcrds.push_back(refloop_[0].ncrd());
	atomcrds.push_back(refloop_[0].cacrd());
	atomcrds.push_back(refloop_[0].ccrd());
	atomcrds.push_back(refloop_[posi].ncrd());
	atomcrds.push_back(refloop_[posi].cacrd());
	atomcrds.push_back(refloop_[posi].ccrd());
	atomcrds.push_back(refloop_[looplength_ - 1].ncrd());
	atomcrds.push_back(refloop_[looplength_ - 1].cacrd());
	atomcrds.push_back(refloop_[looplength_ - 1].ccrd());
	NSPloopclosure::LoopSolver::getSolutions(atomcrds, fixcrds_, solutions_);
	return solutions_.size();
}
