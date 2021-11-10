/*
 * backboneloop.cpp
 *
 *  Created on: 2016年12月7日
 *      Author: hyliu
 */

#include "backbone/backboneloop.h"
#include "dstl/randomengine.h"
#include "backbone/loopinteractions.h"
#include "dstl/sortindex.h"

using namespace NSPproteinrep;
/*
 BackBoneLoop::BackBoneLoop(const BackBoneSite &first, const BackBoneSite &last,
 unsigned int length):first_(first), last_(last) {
 length_ = length + 2;
 chaintopo_ = makeMainChainHeavyTopo(length_);
 crd_ = std::shared_ptr < ChainTreeCrd
 > (new ChainTreeCrd(chaintopo_.get()));
 initConf(first,last);
 isgly_.resize(length_, false);
 ispro_.resize(length_, false);
 }
 */
BackBoneLoop::BackBoneLoop(std::vector<std::vector<BackBoneSite>> *context,
		unsigned int headsegment, unsigned int tailsegment, unsigned int length,
		const std::vector<BackBoneSite> & loopsegment) {
	topsolutions_.initN(nsolns_save);
	contextsegments_ = context;
	headsegment_ = headsegment;
	tailsegment_ = tailsegment;
	length_ = length + 2;
	chaintopo_ = makeMainChainHeavyTopo(length_);
	crd_ = std::shared_ptr < ChainTreeCrd
			> (new ChainTreeCrd(chaintopo_.get()));
	first_ = context->at(headsegment_).back();
	last_ = context->at(tailsegment_)[0];
	BackBoneSite & before =
			context->at(headsegment_)[context->at(headsegment_).size() - 2];
	if (before.omiga() > -90.0 && before.omiga() < 90.0)
		firstcispro_ = true;
	initConf(first_, last_);
	ends_.r_c0 = before.getcrd(BackBoneSite::CCRD);
	ends_.r_ca3 = last_.getcrd(BackBoneSite::CACRD);
	ends_.r_c3 = last_.getcrd(BackBoneSite::CCRD);
	ends_.r_o3 = last_.getcrd(BackBoneSite::OCRD);
	ends_.r_n4 = context->at(tailsegment_)[1].getcrd(BackBoneSite::NCRD);
	isgly_.resize(length_, false);
	ispro_.resize(length_, false);
	if (!loopsegment.empty()) {
		assert(loopsegment.size() == length);

		if ((*context)[headsegment_].back().resname == "GLY")
			isgly_[0] = true;
		else if ((*context)[headsegment_].back().resname == "PRO")
			ispro_[0] = true;
		for (unsigned int i = 0; i < loopsegment.size(); ++i) {
			isgly_[i + 1] = loopsegment[i].resname == "GLY";
			ispro_[i + 1] = loopsegment[i].resname == "PRO";
		}
		if ((*context)[tailsegment_][0].resname == "GLY")
			isgly_[length_ - 1] = true;
		else if ((*context)[tailsegment_][0].resname == "PRO")
			ispro_[length_ - 1] = true;
		chooseglypro_ = false;
	}
}
void BackBoneLoop::initConf(const BackBoneSite &first,
		const BackBoneSite &last) {
	crd_->resetMaps();
	crd_->initWithIdealIntCrd();
	double deg = 3.14159265358979323846 / 180.0;
	double omiga = 3.14159265358979323846;
	if (firstcispro_)
		omiga = 0;
	crd_->setPhi(0, first.phi() * deg);
	crd_->setPsi(0, first.psi() * deg);
	crd_->setOmega(0, omiga);
//	crd_->setPhi(length_ - 1, last.phi() * deg);
//	crd_->setPsi(length_ - 1, last.psi() * deg);
	AtomKey n_k = chaintopo_->genKey(0,
			BackBoneData::atomnames[BackBoneData::N]);
	AtomKey ca_k = chaintopo_->genKey(0,
			BackBoneData::atomnames[BackBoneData::CA]);
	AtomKey c_k = chaintopo_->genKey(0,
			BackBoneData::atomnames[BackBoneData::C]);
	AtomKey o_k = chaintopo_->genKey(0,
			BackBoneData::atomnames[BackBoneData::O]);
	NSPgeometry::XYZ ncrd = first.getcrd(BackBoneSite::NCRD);
	NSPgeometry::XYZ ccrd = first.getcrd(BackBoneSite::CCRD);
	NSPgeometry::XYZ cacrd = first.getcrd(BackBoneSite::CACRD);
	NSPgeometry::XYZ ocrd = first.getcrd(BackBoneSite::OCRD);
	std::map<AtomKey, NSPgeometry::XYZ> crd;
	crd.insert(std::make_pair(n_k, ncrd));
	crd.insert(std::make_pair(ca_k, cacrd));
	crd.insert(std::make_pair(c_k, ccrd));
	crd.insert(std::make_pair(o_k, ocrd));
	crd_->copyCrdMap(crd);
	ends_.r_n1 = ncrd;
	ends_.r_ca1 = cacrd;
	ends_.r_c1 = ccrd;
	ends_.r_o1 = ocrd;

}
void BackBoneLoop::chooseglypro() {
	NSPdstl::RandomEngine<> &rneg = NSPdstl::RandomEngine<>::getinstance();
	const double PGLY { 0.024 };
	const double PPRO { 0.049 };
	rneg.setrealrng(0, 1);
	for (int i = 0; i < length_; ++i) {
		if (i == 0 && firstcispro_) {
			ispro_[i] = true;
			continue;
		}
		isgly_[i] = false;
		ispro_[i] = false;
		crd_->setOmega(i, 3.14159265);
		double rn = rneg.realrng()();
		if (rn < PGLY) {
			isgly_[i] = true;
		} else if (rn < PGLY + PPRO) {
			ispro_[i] = true;
			if (i == 0)
				continue;
			double q = rneg.realrng()();
			if (q > 0.90)
				crd_->setOmega(i, 0);
		}
	}
}
const NSPpdbstatistics::PhiPsiDistr * BackBoneLoop::phipsidistr(
		unsigned int posi) const {
	const NSPpdbstatistics::PhiPsiDistr *distr;
	distr = &(NSPpdbstatistics::PhiPsiDistr::coildistr());
	if (isgly_[posi]) {
		distr = &(NSPpdbstatistics::PhiPsiDistr::glydistr());
	} else if (ispro_[posi]) {
		double w = crd_->getOmega(posi);
		if (w > -0.1 && w < 0.1)
			distr = &(NSPpdbstatistics::PhiPsiDistr::cisprodistr());
		else
			distr = &(NSPpdbstatistics::PhiPsiDistr::transprodistr());

	} else if (posi < length_ - 1) {
		if (ispro_[posi + 1])
			distr = &(NSPpdbstatistics::PhiPsiDistr::preprodistr());
	}
	return distr;
}

void BackBoneLoop::randomConf(bool choosegp) {
	double rad = 180.0 / 3.14159265;
	if (choosegp)
		chooseglypro();
	NSPdstl::RandomEngine<> &rneg = NSPdstl::RandomEngine<>::getinstance();
	for (int posi = 0; posi < length_; ++posi) {
		const NSPpdbstatistics::PhiPsiDistr *distr = phipsidistr(posi);
		double phi, psi;
		if (posi == 0) {
			phi = crd_->getPhi(posi) * rad;
			distr->randompsi(rneg.realrng(), phi, &psi);
		} else {
			distr->randomphipsi(rneg.realrng(), &phi, &psi);
		}
		crd_->setPhi(posi, phi / rad);
		crd_->setPsi(posi, psi / rad);
	}
	crd_->calcXYZ(true);
}
/*
void BackBoneLoop::partial_randomConf() {
	double rad = 180.0 / 3.14159265;
	NSPdstl::RandomEngine<> &rneg = NSPdstl::RandomEngine<>::getinstance();
	double pchange = 0.1 + rneg.realrng()() * 0.2;
	; // 1/10 to 3/10 positions will be changed
	std::vector<int> tochange;
	for (int posi = 0; posi < length_; ++posi) {
		if (rneg.realrng()() >= pchange)
			continue;
		tochange.push_back(posi);
	}
	if (tochange.empty())
		tochange.push_back(rneg.intrng(0, (int) (length_ - 1))());
	for (int posi : tochange) {
		const NSPpdbstatistics::PhiPsiDistr *distr = phipsidistr(posi);
		double phi, psi;
		if (posi == 0) {
			phi = crd_->getPhi(posi) * rad;
			distr->randompsi(rneg.realrng(), phi, &psi);
		} else {
			distr->randomphipsi(rneg.realrng(), &phi, &psi);
		}
		crd_->setPhi(posi, phi / rad);
		crd_->setPsi(posi, psi / rad);
	}
	crd_->calcXYZ(true);
}
*/
double BackBoneLoop::scoreSolution(const ClosureSolution &sol) const {
	double res = 0.0;
	const NSPpdbstatistics::PhiPsiDistr *distr = phipsidistr(sol.posi1);
	double rad = 180.0 / 3.14159265;
	double phi1 = crd_->getPhi(sol.posi1) * rad;
	res += distr->statisticalenergy(phi1, sol.psi1);
	distr = phipsidistr(sol.posi2);
	res += distr->statisticalenergy(sol.phi2, sol.psi2);
	distr = phipsidistr(sol.posi3);
	double psi3 = crd_->getPsi(sol.posi3);
	res += distr->statisticalenergy(sol.phi3, psi3);
	return res;
}

double BackBoneLoop::torsionscore() const {
	double res = 0.0;
	static const double rad = 180.0 / 3.14159265358979323846;
	for (unsigned int posi = 0; posi < length_; ++posi) {
		const NSPpdbstatistics::PhiPsiDistr *distr = phipsidistr(posi);
		res += distr->statisticalenergy(crd_->getPhi(posi) * rad,
				crd_->getPsi(posi) * rad);
	}
//	std::cout <<"XXXLoop torsion score: " << res <<std::endl;
	return res;
}
int BackBoneLoop::batchFindSolutions(unsigned int ncopies) {
	NSPdstl::RandomEngine<> &rneg = NSPdstl::RandomEngine<>::getinstance();
	int ntry = 0;
	int nsuccess=0;
	while (nsuccess < ncopies && ntry < 50 * ncopies) {
		initConf(first_, last_);
		randomConf(chooseglypro_);
		std::vector<CrdSolution> solns;
		for (unsigned int posi = 1; posi < length_ - 1; ++posi)
			closureSolutions(posi, solns);
		if (!solns.empty())
			nsuccess +=saveTopSolutions(solns);
		++ntry;
	}
	return nsuccess;
}
int BackBoneLoop::saveTopSolutions(std::vector<CrdSolution> &solns) {
	int nsuccess = 0;
	for (auto & s : solns) {
		double score = scoreCrdSolution(s);
		if(score > LoopInteractions::SCORECUT) continue;
		++nsuccess;
		bool save = score < refscore_ && saveacceptable_;
		if (topsolutions_.keep(score) || save) {
			auto loop = std::shared_ptr < std::vector
					< BackBoneSite >> (new std::vector<BackBoneSite>);
			makeAllsites(loop.get(), 1);
			topsolutions_.push(loop, score);
			if (save) {
				acceptablesolutions_.push_back(loop);
				acceptablescores_.push_back(score);
			}
		}
	}
	return nsuccess;
}

bool BackBoneLoop::newLoopConf(std::vector<BackBoneSite> *loop, double *score) {
	NSPdstl::RandomEngine<> &rneg = NSPdstl::RandomEngine<>::getinstance();
	bool res = false;
	unsigned ntry = 0;
	while (!res && ntry++ < 100) {
		initConf(first_, last_);
		randomConf(chooseglypro_);
		std::vector<CrdSolution> solns;
		for (unsigned int posi = 1; posi < length_ - 1; ++posi)
			closureSolutions(posi, solns);
		if (!solns.empty()) {
			*score = chooseApplySolution(solns, 0.001);
			if (*score < 10000.0) {
				res = true;
			}
		}
		if (res)
			makeAllsites(loop, 1);
	}
//	std::cout <<"Ntry: " <<ntry <<std::endl;
	return res;
}

unsigned int BackBoneLoop::getLoops(unsigned int ncopies,
		unsigned int ncandidates,
		typename std::vector<std::vector<BackBoneSite>>::iterator loopiter,
		std::vector<double> *scores) {
	unsigned int m = 0;
	unsigned int n=0;
	unsigned int ntry = 0;
	unsigned int mintry = ncandidates;
	unsigned int maxtry = 200 * mintry;
	double scoremax=0;
	topsolutions_.initN(ncopies);
	while (n<ncopies || (m < ncopies/2 && ntry < maxtry)) {
		batchFindSolutions(ncopies);
		++ntry;
		if (ntry % 50 == 0) {
			std::cout <<"After " << ntry<< " batches, "
					<< n <<" copies of top solutions saved, with scores below " << scoremax
					<<". " << m
					<< " copies of acceptable loops generated." << std::endl;
		}
		m=acceptablesolutions_.size();
		if(saveacceptable_) {
			if( ntry > 300 && m <= 1) break;
		}
	    n=topsolutions_.size();
	    if(n>0) topsolutions_.top(&scoremax);
	}
	auto topn=NSPdstl::topN2vector(topsolutions_);
	if(scores) {
		scores->resize(topn.size(),0.0);
	}
	n=0;
	for(auto iter=topn.begin(); iter != topn.end();++iter){
		*(loopiter++) = *(iter->first);
		(*scores)[n++]=iter->second;
	}
	return n;
}

/*
 unsigned int BackBoneLoop::getLoops(unsigned int ncopies,
 unsigned int ncandidates,
 typename std::vector<std::vector<BackBoneSite>>::iterator loopiter,
 std::vector<double> *scores) {
 NSPdstl::TopN<std::shared_ptr<std::vector<BackBoneSite>>>toploops(ncopies);
 unsigned int m=0;
 unsigned int ntry=0;
 unsigned int mintry=ncopies*ncandidates;
 unsigned int maxtry=100*mintry;
 while(ntry < mintry ||(m <ncopies && ntry <maxtry)) {
 ++ntry;
 if((ntry % 100) == 0) {
 std::cout <<m <<" copies of loop generated after " <<ntry << " tries." <<std::endl;
 }
 auto loop=std::shared_ptr<std::vector<BackBoneSite>> (new std::vector<BackBoneSite>);
 double score;
 if(!newLoopConf(loop.get(),&score)) continue;
 if(score > refscore_) continue;
 toploops.push(loop,score);
 m=toploops.size();
 }
 if(m>ncopies) m=ncopies;
 double score;
 if(scores) scores->resize(m,0.0);
 for(unsigned int n=0;n<m; ++n) {
 *(loopiter++)=*(toploops.top(&score));
 if(scores) (*scores)[n]=score;
 //			std::cout <<"Loop score: " << -score <<std::endl;
 toploops.pop();
 }
 return m;
 }
*/

double BackBoneLoop::copyNScoreLoop(std::vector<BackBoneSite>::const_iterator &iter) {
	crd_->copyCrdFromSites(iter,firstcispro_);
	if(chooseglypro_){
		for(int i=0; i<length_;++i) {
			isgly_[i]=false;
			ispro_[i]=false;
			auto bs=iter+i;
			if(bs->resname=="GLY") isgly_[i]=true;
			else if(bs->resname=="PRO") ispro_[i]=true;
		}
	}
	return loopScore();
}
double BackBoneLoop::loopScore() {
	double ene = 0.0;
	if (contextsegments_) {
		std::vector<BackBoneSite> loop;
		makeAllsites(&loop, 0);
		ene = loop_context_interaction(contextsegments_, headsegment_,
				tailsegment_, loop);
	}
	double ts = torsionscore();
	return ene + ts;
}
double BackBoneLoop::scoreCrdSolution(const CrdSolution &s) {
	crd_->copyCrdMap(s.crdmap); // internal torsional anlges will be updated
	crd_->setPhi(s.posi1, s.phi1);  //phi1 cannot be calculated from coordinates
	crd_->setPsi(s.posi3, s.psi3);
	return loopScore();
}

double BackBoneLoop::chooseApplySolution(std::vector<CrdSolution> &solns,
		double temperature) {
	std::vector<double> fact;
	std::vector<double> score;
	double facttot;
	for (auto & s : solns) {
		double sc = scoreCrdSolution(s);
		score.push_back(sc);
		double bfac = exp(-sc / temperature);
		facttot += bfac;
		fact.push_back(bfac);
	}
	NSPdstl::RandomEngine<> &rneg = NSPdstl::RandomEngine<>::getinstance();
	double rn = rneg.realrng()();
	double facsum = 0.0;
	double sc = score.back();
	CrdSolution *s = &(solns.back());
	for (unsigned int i = 0; i < fact.size() - 1; ++i) {
		facsum += fact[i] / facttot;
		if (rn < facsum) {
			s = &(solns.at(i));
			sc = score[i];
			break;
		}
	}
	crd_->copyCrdMap(s->crdmap);
	crd_->setPhi(s->posi1, s->phi1);
	crd_->setPsi(s->posi3, s->psi3);
	return sc;
//	crd_->calcXYZ();
}

void BackBoneLoop::closureSolutions(unsigned int posi,
		std::vector<CrdSolution> &solns) const {
	assert(posi > 0 && posi < length_ - 1);
	std::vector<NSPgeometry::XYZ> fixcrds = ends_.fixedCrds();
	std::vector<NSPgeometry::XYZ> atomcrds;
	atomcrds.push_back(
			crd_->getCrd(0, BackBoneData::atomnames[BackBoneData::N]));
	atomcrds.push_back(
			crd_->getCrd(0, BackBoneData::atomnames[BackBoneData::CA]));
	atomcrds.push_back(
			crd_->getCrd(0, BackBoneData::atomnames[BackBoneData::C]));
	atomcrds.push_back(
			crd_->getCrd(posi, BackBoneData::atomnames[BackBoneData::N]));
	atomcrds.push_back(
			crd_->getCrd(posi, BackBoneData::atomnames[BackBoneData::CA]));
	atomcrds.push_back(
			crd_->getCrd(posi, BackBoneData::atomnames[BackBoneData::C]));
	atomcrds.push_back(
			crd_->getCrd(length_ - 1,
					BackBoneData::atomnames[BackBoneData::N]));
	atomcrds.push_back(
			crd_->getCrd(length_ - 1,
					BackBoneData::atomnames[BackBoneData::CA]));
	atomcrds.push_back(
			crd_->getCrd(length_ - 1,
					BackBoneData::atomnames[BackBoneData::C]));
	std::vector<NSPloopclosure::LoopSolution> solutions;
	NSPloopclosure::LoopSolver::getSolutions(atomcrds, fixcrds, solutions);
	for (auto & sol : solutions) {
		solns.push_back(
				CrdSolution(ends_, crd_.get(), 0, posi, length_ - 1, sol));
	}
}
void BackBoneLoop::closureSolutions(unsigned int p1, unsigned int p2,
		unsigned int p3, std::vector<ClosureSolution> &solns) const {
	std::vector<NSPgeometry::XYZ> fixcrds;
	std::vector<NSPgeometry::XYZ> atomcrds;
	atomcrds.push_back(
			crd_->getCrd(p1, BackBoneData::atomnames[BackBoneData::N]));
	atomcrds.push_back(
			crd_->getCrd(p1, BackBoneData::atomnames[BackBoneData::CA]));
	atomcrds.push_back(
			crd_->getCrd(p1, BackBoneData::atomnames[BackBoneData::C]));
	atomcrds.push_back(
			crd_->getCrd(p2, BackBoneData::atomnames[BackBoneData::N]));
	atomcrds.push_back(
			crd_->getCrd(p2, BackBoneData::atomnames[BackBoneData::CA]));
	atomcrds.push_back(
			crd_->getCrd(p2, BackBoneData::atomnames[BackBoneData::C]));
	atomcrds.push_back(
			crd_->getCrd(p3, BackBoneData::atomnames[BackBoneData::N]));
	atomcrds.push_back(
			crd_->getCrd(p3, BackBoneData::atomnames[BackBoneData::CA]));
	atomcrds.push_back(
			crd_->getCrd(p3, BackBoneData::atomnames[BackBoneData::C]));
	fixcrds.push_back(atomcrds[0]);
	fixcrds.push_back(atomcrds[1]);
	fixcrds.push_back(atomcrds[7]);
	fixcrds.push_back(atomcrds[8]);
	std::vector<NSPloopclosure::LoopSolution> solutions;
	NSPloopclosure::LoopSolver::getSolutions(atomcrds, fixcrds, solutions);
	for (auto & sol : solutions) {
		solns.push_back(ClosureSolution(ends_, crd_.get(), p1, p2, p3, sol));
	}
}

/*struct CrdSolution {
 std::map<typename ChainTreeTopology::AtomKey,NSPgeometry::XYZ> crdmap;*/
CrdSolution::CrdSolution(const FixEnds & ends, ChainTreeCrd *crd,
		unsigned int p1, unsigned int p2, unsigned int p3,
		const NSPloopclosure::LoopSolution &sol) :
		posi1(p1), posi2(p2), posi3(p3) {
	const NSPgeometry::RigidTransform & rt = sol.rt1;
	typedef typename ChainTreeTopology::AtomKey AtomKey;
	for (unsigned int p = p1; p <= p2; ++p) {
		AtomKey nk = crd->topo()->genKey(p,
				BackBoneData::atomnames[BackBoneData::N]);
		AtomKey cak = crd->topo()->genKey(p,
				BackBoneData::atomnames[BackBoneData::CA]);
		AtomKey ck = crd->topo()->genKey(p,
				BackBoneData::atomnames[BackBoneData::C]);
		AtomKey ok = crd->topo()->genKey(p,
				BackBoneData::atomnames[BackBoneData::O]);
		if (p == 0) {
			crdmap.insert(std::make_pair(nk, ends.r_n1));
			crdmap.insert(std::make_pair(cak, ends.r_ca1));
		} else if (p > p1) {
			crdmap.insert(
					std::make_pair(nk,
							rt.applytoCopy(
									crd->getCrd(p,
											BackBoneData::atomnames[BackBoneData::N]))));
			crdmap.insert(
					std::make_pair(cak,
							rt.applytoCopy(
									crd->getCrd(p,
											BackBoneData::atomnames[BackBoneData::CA]))));
		}
		if (p != p2) {
			crdmap.insert(
					std::make_pair(ck,
							rt.applytoCopy(
									crd->getCrd(p,
											BackBoneData::atomnames[BackBoneData::C]))));
			crdmap.insert(
					std::make_pair(ok,
							rt.applytoCopy(
									crd->getCrd(p,
											BackBoneData::atomnames[BackBoneData::O]))));
		}
	}
	NSPgeometry::XYZ c1 = rt.applytoCopy(
			crd->getCrd(p1, BackBoneData::atomnames[BackBoneData::C]));
	phi1 = NSPgeometry::torsion(ends.r_c0, ends.r_n1, ends.r_ca1, c1);
	const NSPgeometry::RigidTransform & rt2 = sol.rt2;
	for (unsigned int p = p2; p <= p3; ++p) {
		AtomKey nk = crd->topo()->genKey(p,
				BackBoneData::atomnames[BackBoneData::N]);
		AtomKey cak = crd->topo()->genKey(p,
				BackBoneData::atomnames[BackBoneData::CA]);
		AtomKey ck = crd->topo()->genKey(p,
				BackBoneData::atomnames[BackBoneData::C]);
		AtomKey ok = crd->topo()->genKey(p,
				BackBoneData::atomnames[BackBoneData::O]);
		if (p > p2) {
			crdmap.insert(
					std::make_pair(nk,
							rt2.applytoCopy(
									crd->getCrd(p,
											BackBoneData::atomnames[BackBoneData::N]))));
		}
		if (p != p3) {
			crdmap.insert(
					std::make_pair(cak,
							rt2.applytoCopy(
									crd->getCrd(p,
											BackBoneData::atomnames[BackBoneData::CA]))));
			crdmap.insert(
					std::make_pair(ck,
							rt2.applytoCopy(
									crd->getCrd(p,
											BackBoneData::atomnames[BackBoneData::C]))));
			crdmap.insert(
					std::make_pair(ok,
							rt2.applytoCopy(
									crd->getCrd(p,
											BackBoneData::atomnames[BackBoneData::O]))));
		} else if (p3 == crd->length() - 1) {
			crdmap.insert(std::make_pair(cak, ends.r_ca3));
			crdmap.insert(std::make_pair(ck, ends.r_c3));
			crdmap.insert(std::make_pair(ok, ends.r_o3));
		}
	}
	NSPgeometry::XYZ n3 = rt2.applytoCopy(
			crd->getCrd(p3, BackBoneData::atomnames[BackBoneData::N]));
	psi3 = NSPgeometry::torsion(n3, ends.r_ca3, ends.r_c3, ends.r_n4);
}
;
ClosureSolution::ClosureSolution(const FixEnds & ends, ChainTreeCrd *crd,
		unsigned int p1, unsigned int p2, unsigned int p3,
		const NSPloopclosure::LoopSolution &sol) :
		posi1(p1), posi2(p2), posi3(p3) {
	NSPgeometry::XYZ c1 = sol.rt1.applytoCopy(
			crd->getCrd(p1, BackBoneData::atomnames[BackBoneData::C]));
	NSPgeometry::XYZ n2 = sol.rt1.applytoCopy(
			crd->getCrd(p1 + 1, BackBoneData::atomnames[BackBoneData::N]));
	phi1 = NSPgeometry::torsion(ends.r_c0, ends.r_n1, ends.r_ca1, c1);
	psi1 = NSPgeometry::torsion(ends.r_n1, ends.r_ca1, c1, n2);
	c1 = sol.rt1.applytoCopy(
			crd->getCrd(p2 - 1, BackBoneData::atomnames[BackBoneData::C]));
	n2 = sol.rt1.applytoCopy(
			crd->getCrd(p2, BackBoneData::atomnames[BackBoneData::N]));
	NSPgeometry::XYZ ca2 = sol.rt1.applytoCopy(
			crd->getCrd(p2, BackBoneData::atomnames[BackBoneData::CA]));
	NSPgeometry::XYZ c2 = sol.rt2.applytoCopy(
			crd->getCrd(p2, BackBoneData::atomnames[BackBoneData::C]));
	NSPgeometry::XYZ n3 = sol.rt2.applytoCopy(
			crd->getCrd(p2 + 1, BackBoneData::atomnames[BackBoneData::N]));
	std::cout << "NCRDS: " << std::endl;
	std::cout << "n2: " << n2.toString() << std::endl;
	std::cout << "ca2: " << ca2.toString() << std::endl;
	std::cout << "c2: " << c2.toString() << std::endl;
	phi2 = torsion(c1, n2, ca2, c2);
	psi2 = torsion(n2, ca2, c2, n3);
	std::cout << "phi2, psi2 " << phi2 * 180.0 / 3.14159265 << " "
			<< psi2 * 180.0 / 3.14159265 << std::endl;
	std::cout << "angles " << angle(n2, ca2, c2) * 180 / 3.14159265 << " "
			<< angle(ca2, c2, n3) * 180 / 3.14159265 << std::endl;
	c2 = sol.rt2.applytoCopy(
			crd->getCrd(p3 - 1, BackBoneData::atomnames[BackBoneData::C]));
	n3 = sol.rt2.applytoCopy(
			crd->getCrd(p3, BackBoneData::atomnames[BackBoneData::N]));
	phi3 = torsion(c2, n3, ends.r_ca3, ends.r_c3);
	psi3 = torsion(n3, ends.r_ca3, ends.r_c3, ends.r_n4);
}

BackBoneSite BackBoneLoop::makeSite(unsigned int posi, unsigned int shift,
		std::string pdbid, char chainid, std::string resname) const {
	BackBoneSite bs;
	bs.pdbid = pdbid;
	bs.chainid = chainid;
	bs.resid = posi + shift;
	bs.resseq = posi + shift;
	bs.sscode = 'C';
	if (isgly_[posi])
		resname = "GLY";
	else if (ispro_[posi])
		resname = "PRO";
	bs.resname = resname;
	std::vector<NSPgeometry::XYZ> sitecrd;
	crd_->getBackBoneCrd(posi, sitecrd);
	bs.changecrd(sitecrd);
	double rad = 180.0 / 3.14159265358979323846;
	double t;
//	if(posi != 0) {
	t = crd_->getPhi(posi) * rad;
	while (t > 180.0)
		t -= 360.0;
	while (t < -180.0)
		t += 360.0;
	bs.data_[BackBoneSite::PHI] = t;
//	} else {
//		BackBoneSite &ps=(*contextsegments_)[headsegment_][(*contextsegments_)[headsegment_].size()-2];
//		bs.phi(ps);
//	}
	if (posi != length_ - 1) {
		t = crd_->getPsi(posi) * rad;
		while (t > 180.0)
			t -= 360.0;
		while (t < -180.0)
			t += 360.0;
		bs.data_[BackBoneSite::PSI] = t;
		t = crd_->getOmega(posi + 1) * rad;
		while (t > 180.0)
			t -= 360.0;
		while (t < -180.0)
			t += 360.0;
		bs.data_[BackBoneSite::OMIGA] = t;
	} else {
		BackBoneSite &ns = (*contextsegments_)[tailsegment_][1];
//			std::cout <<"psi-compare: " <<bs.psi();
//		bs.psi(ns);
		t = crd_->getPsi(posi) * rad;
		while (t > 180.0)
			t -= 360.0;
		while (t < -180.0)
			t += 360.0;
		bs.data_[BackBoneSite::PSI] = t;
		bs.omiga(ns);
	}
	return bs;
}
