/*
 * closingloop.cpp
 *
 *  Created on: 2017年4月23日
 *      Author: hyliu
 */

#include "backbone/mainchain.h"
#include "backbone/closingloop.h"
#include "pdbstatistics/phipsidistr.h"
#include "dstl/randomengine.h"
#include <algorithm>
using namespace NSPproteinrep;

ClosingLoop::ClosingLoop(const MainChain &host, int startposi, int endposi,
		int modestype, int pertposi) {
	assert(startposi >= 1);
	assert(endposi < host.size());
	int length = endposi - startposi;
	setmodes(modestype);
	if (confmode_ == MUTATE) {
		if (pertposi >= 0)
			pertposi_ = pertposi - startposi;
	}
	for (int p2 = 1; p2 < length - 1; ++p2) {
		if (host.localrigid(startposi+p2))
			continue;
		allowedp2_.push_back(p2);
	}
	if (host.localrigid(startposi) || host.localrigid(endposi - 1)) {
		allowedp2_.clear();
	}
	tmplt_.resize(length, BackBoneSite());
	host.copysegment(startposi, length, &tmplt_);
	fixcrds_.push_back(tmplt_[0].ncrd());
	fixcrds_.push_back(tmplt_[0].cacrd());
	fixcrds_.push_back(tmplt_.back().cacrd());
	fixcrds_.push_back(tmplt_.back().ccrd());
	rc0_ = host[startposi - 1].ccrd();
	rn4_ = host[endposi].ncrd();
	lastomiga_=tmplt_.back().omiga();
	if (endposi < host.size())
		nextpro_ = host[endposi].resname == "PRO";

}
int ClosingLoop::solve() {
	if (allowedp2_.empty())
		return 0;
	initseqconf();
//	for (auto i = refsites_.begin(); i != refsites_.end(); ++i) {
//		std::cout << i->phi() << "\t" << i->psi() << "\t";
//	}
//	std::cout << std::endl;
	clearsolutions();
// just for test
	int nsol_old = solutions_.size();
	for (auto p2 : allowedp2_) {
//		if(confmode_== MUTATE && p2==selectedpertposi_) continue;
		int nsol_new = addsolutions(p2);
		for (int i = nsol_old; i < nsol_new; ++i) {
			p2ofsolutions_.push_back(p2);
		}
		nsol_old=nsol_new;
	}
	return solutions_.size();
}
void ClosingLoop::getsitessolution(int i, std::vector<BackBoneSite> *result) {
	NSPloopclosure::LoopSolution &sol = solutions_.at(i);
	int p2 = p2ofsolutions_.at(i);
	const NSPgeometry::RigidTransform & rt1 = sol.rt1;
	const NSPgeometry::RigidTransform & rt2 = sol.rt2;
	int length = refsites_.size();
	result->resize(length, BackBoneSite());
	std::copy(refsites_.begin(), refsites_.end(), result->begin());

	for (int p = 0; p < length; ++p) {
		std::vector<NSPgeometry::XYZ> crd;
		if (p == 0) {
			crd.push_back(refsites_[p].ncrd());
			crd.push_back(refsites_[p].cacrd());
			crd.push_back(rt1.applytoCopy(refsites_[p].ccrd()));
			crd.push_back(rt1.applytoCopy(refsites_[p].ocrd()));
		} else if (p < p2) {
			crd.push_back(rt1.applytoCopy(refsites_[p].ncrd()));
			crd.push_back(rt1.applytoCopy(refsites_[p].cacrd()));
			crd.push_back(rt1.applytoCopy(refsites_[p].ccrd()));
			crd.push_back(rt1.applytoCopy(refsites_[p].ocrd()));
		} else if (p == p2) {
			double rca_c1=NSPgeometry::distance(refsites_[p].cacrd(),
					refsites_[p].ccrd());
			crd.push_back(rt1.applytoCopy(refsites_[p].ncrd()));
			crd.push_back(rt1.applytoCopy(refsites_[p].cacrd()));
//
			crd.push_back(rt2.applytoCopy(refsites_[p].ccrd()));
			crd.push_back(rt2.applytoCopy(refsites_[p].ocrd()));
			double rca_c2=NSPgeometry::distance(rt1.applytoCopy(refsites_[p].cacrd()),
					rt2.applytoCopy(refsites_[p].ccrd()));
			if(rca_c2>1.6) {
				std::cout<<"rca_c at p2: " <<rca_c1 <<"\t" <<rca_c2<<std::endl;
			}
		} else if (p != length - 1) {
			crd.push_back(rt2.applytoCopy(refsites_[p].ncrd()));
			crd.push_back(rt2.applytoCopy(refsites_[p].cacrd()));
			crd.push_back(rt2.applytoCopy(refsites_[p].ccrd()));
			crd.push_back(rt2.applytoCopy(refsites_[p].ocrd()));
		} else {
			crd.push_back(rt2.applytoCopy(refsites_[p].ncrd()));
			crd.push_back(fixcrds_[2]);
			crd.push_back(fixcrds_[3]);
			NSPgeometry::XYZ ocrd=NSPgeometry::InternaltoXYZ(fixcrds_[3],
					fixcrds_[2],rn4_,1.23,2.164,3.14159265);
			crd.push_back(ocrd);
		}
		result->at(p).changecrd(crd.begin());
		result->at(p).sscode='C';
	}
	double phi1, psi1, phi2, psi2, phi3, psi3;
	double rad = 180.0 / 3.14159265;
	BackBoneSite &s1 = result->at(0);
	phi1 = NSPgeometry::torsion(rc0_, s1.ncrd(), s1.cacrd(), s1.ccrd());
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
	psi3 = NSPgeometry::torsion(s3.ncrd(), s3.cacrd(), s3.ccrd(), rn4_);
	s3.data_[BackBoneSite::PHI] = phi3 * rad;
	s3.data_[BackBoneSite::PSI] = psi3 * rad;
	s3.data_[BackBoneSite::OMIGA]=lastomiga_;
/*
	for(int m=1;m<result->size()-1;++m){
		BackBoneSite &bs=result->at(m);
		BackBoneSite &ps=result->at(m-1);
		BackBoneSite &ns=result->at(m+1);
		double phi=bs.phi();
		double psi=bs.psi();
		double omiga=ps.omiga();
		double rphi=bs.phi(ps);
		double rpsi=bs.psi(ns);
		double romiga=ps.omiga(bs);
		double diff_p=phi-rphi;
		double diff_psi=psi-rpsi;
		double diff_o=omiga-romiga;
		auto shift=[](double t){if(t>180.0) t-=360.0; if(t<-180.0) t+=360.0;return t;};
		auto large=[](double d,double d2=1.e-5){return d>d2 || d<-d2;};
		if(large(shift(diff_p)) || large(shift(diff_psi))||large(shift(diff_o),10.0)){
			std::cout <<"Inconsistent loop site "<<diff_p <<"\t"<<diff_psi <<"\t"<<diff_o<<std::endl;
			exit(1);
		}
		std::cout <<m+1<<"\t"<<phi <<"\t" <<psi <<"\t" <<omiga<<std::endl;
	}
	*/
}

void ClosingLoop::initseqconf() {
	selectedpertposi_ = pertposi_;
	if ((seqmode_ == MUTATE || confmode_==MUTATE) && pertposi_ <= 0)
		selectedpertposi_ = NSPdstl::RandomEngine<>::getinstance().intrng(1,
				tmplt_.size() - 2)();
	initseq();
	initconf();
}
void ClosingLoop::initseq() {
	int length = tmplt_.size();
	seq_.resize(length, "ALA");
	if (seqmode_ == NEW) {
		for (int i = 0; i < length; i++)
			seq_[i] = chooseresiduetype();
	} else {
		for (int i = 0; i < length; ++i) {
			if (tmplt_.missing(i)
					|| (seqmode_ == MUTATE && i == selectedpertposi_))
				seq_[i] = chooseresiduetype();
			else
				seq_[i] = tmplt_[i].resname;
		}
	}
}
void ClosingLoop::initconf() {
	int length = tmplt_.size();
	refsites_.resize(length, BackBoneSite());
	refsites_[0] = tmplt_[0];
	for (int i = 1; i < length; ++i) {
		double phi = tmplt_[i].phi();
		double psi = tmplt_[i].psi();
		double omiga= refsites_[i-1].omiga();
		if(!tmplt_.missing(i-1))  omiga=tmplt_[i-1].omiga();
		bool prevcis = omiga>-90.0 &&omiga <90.0;
		std::string next = "ALA";
		if (i < length - 1)
				next = seq_[i + 1];
		else if (nextpro_)
				next = "PRO";
		if (confmode_ == NEW || tmplt_.missing(i)) {
			choosephipsi(seq_[i], next, &phi, &psi, &prevcis);
		} else if(i == selectedpertposi_) {
			if(NSPdstl::RandomEngine<>::getinstance().realrng(0,1.0)() <0.01)
				choosephipsi(seq_[i], next, &phi, &psi, &prevcis);
			else {
				double r1=NSPdstl::RandomEngine<>::getinstance().realrng(0,1.0)();
				phi +=  -5.0+r1*10.0;
				r1=NSPdstl::RandomEngine<>::getinstance().realrng(0,1.0)();
				psi += -5.0+r1*10.0;
			}
		}
		refsites_[i].resname = seq_[i];
		genbackbonesite(&refsites_[i - 1], prevcis, phi, psi, &refsites_[i]);
	}
	refsites_[length-1].data_[BackBoneSite::OMIGA]=lastomiga_;
}
void ClosingLoop::setmodes(int modestype) {
	switch (modestype) {
	case (NEWSEQ_NEWCONF): {
		seqmode_ = NEW;
		confmode_ = NEW;
		break;
	}
	case (NEWSEQ_MUTATECONF): {
		seqmode_ = NEW;
		confmode_ = MUTATE;
		break;
	}
	case (KEEPSEQ_NEWCONF): {
		seqmode_ = KEEP;
		confmode_ = MUTATE;
		break;
	}
	case (KEEPSEQ_MUTATECONF): {
		seqmode_ = KEEP;
		confmode_ = MUTATE;
		break;
	}
	case (MUTATESEQ_CONF): {
		seqmode_ = MUTATE;
		confmode_ = MUTATE;
		break;
	}
	default: {
		;
	}
	}
}
std::string NSPproteinrep::chooseresiduetype() {
	const double PGLY { 0.122 };
	const double PPRO { 0.086 };
	NSPdstl::RandomEngine<> &rneg = NSPdstl::RandomEngine<>::getinstance();
	std::string res = "ALA";
	rneg.setrealrng(0, 1);
	double rn = rneg.realrng()();
	if (rn < PGLY) {
		res = "GLY";
	} else if (rn < PGLY + PPRO) {
		res = "PRO";
	}
	return res;
}
void NSPproteinrep::choosephipsi(std::string restype,
		std::string nextresiduetype, double *phi, double *psi, bool *cispro) {
	NSPdstl::RandomEngine<> &rneg = NSPdstl::RandomEngine<>::getinstance();
	const double PCIS { 0.1 };

	const NSPpdbstatistics::PhiPsiDistr *distr =
			&(NSPpdbstatistics::PhiPsiDistr::coildistr());
	if (cispro) *cispro=false;
	if (restype == "GLY") {
		distr = &(NSPpdbstatistics::PhiPsiDistr::glydistr());
	} else if (restype == "PRO") {
		distr = &(NSPpdbstatistics::PhiPsiDistr::transprodistr());
		if (cispro) {
			if (rneg.realrng()() < PCIS) {
				*cispro = true;
				distr = &(NSPpdbstatistics::PhiPsiDistr::cisprodistr());
			}
		}
	} else {
		if (nextresiduetype == "PRO") {
			distr = &(NSPpdbstatistics::PhiPsiDistr::preprodistr());
		}
	}
	distr->randomphipsi(rneg.realrng(), phi, psi);
}
int ClosingLoop::addsolutions(int posi) {
	int length = refsites_.size();
	assert(posi > 0 && posi < length - 1);
	std::vector<NSPgeometry::XYZ> atomcrds;
	atomcrds.push_back(refsites_[0].ncrd());
	atomcrds.push_back(refsites_[0].cacrd());
	atomcrds.push_back(refsites_[0].ccrd());
	atomcrds.push_back(refsites_[posi].ncrd());
	atomcrds.push_back(refsites_[posi].cacrd());
	atomcrds.push_back(refsites_[posi].ccrd());
	atomcrds.push_back(refsites_[length - 1].ncrd());
	atomcrds.push_back(refsites_[length - 1].cacrd());
	atomcrds.push_back(refsites_[length - 1].ccrd());
	NSPloopclosure::LoopSolver::getSolutions(atomcrds, fixcrds_, solutions_);
	return solutions_.size();
}

