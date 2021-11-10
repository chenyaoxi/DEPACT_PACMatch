/*
 * backboneenergy.cpp
 *
 *  Created on: 2017年8月2日
 *      Author: hyliu
 */

#include "backbone/backboneenergy.h"
#include "pdbstatistics/phipsidistr.h"
#include "pdbstatistics/proteinblock.h"
#include "pdbstatistics/pbtetrabase.h"
#include "pdbstatistics/nn_pepscorer.h"
#include "dstl/stlutil.h"
#include "backbone/chainpack.h"
#include <algorithm>
#include <map>
using namespace NSPproteinrep;
using namespace NSPdataio;
int EnergyComponents::NUM_ECOMP = { 6 };
//bool print{false};
int NSPproteinrep::initenergycontrols( const std::string &ecid,
		const std::vector<std::string> &controlines) {
	std::map<std::string, double> default_doublepars { { "PhiPsiWeight", 1.0 },
			{ "PBLocalWeight", 1.0 }, { "ClashWeight", 1.0 }, {
					"PBPackingWeight", 1.0 }, { "ClashEnergy", 10000.0 } };
	std::map<std::string, std::string> default_stringpars { { "RefPBSeq", "" },
			{ "PositionMask", "" } };
	std::map<std::string, int> default_intpars { { "PackingSeparation", 6 }, {
			"SSasPBType", 0 },{"IgnoreResName",0} };
	ParameterSet &controlpar = EnergyControls::getparameterset(ecid);
	controlpar.InitDoubleKeys(NSPdstl::getkeyvec(default_doublepars));
	controlpar.InitIntKeys(NSPdstl::getkeyvec(default_intpars));
	controlpar.InitStringKeys(NSPdstl::getkeyvec(default_stringpars));
	controlpar.initialized=true;
	int nread = adjustenergycontrols(ecid,controlines);
	for (auto &kv : default_doublepars) {
		if (!controlpar.keydefined(kv.first)) {
			controlpar.DoublePar.insert(kv.first, kv.second);
		}
	}
	for (auto &kv : default_intpars) {
		if (!controlpar.keydefined(kv.first)) {
			controlpar.IntPar.insert(kv.first, kv.second);
		}
	}
	for (auto &kv : default_stringpars) {
		if (!controlpar.keydefined(kv.first)) {
			controlpar.StringPar.insert(kv.first, kv.second);
		}
	}

	return nread;
}
int NSPproteinrep::adjustenergycontrols(const std::string &ecid,
		const std::vector<std::string> &controlines) {
	ParameterSet & controlpar = EnergyControls::getparameterset(ecid);
	if (!controlpar.initialized)
		initenergycontrols(ecid);
	int nsuccess = 0;
	for (auto &line : controlines)
		if (controlpar.readline(line))
			nsuccess++;
	return nsuccess;
}
void ChainEnergyControl::setpositionmask(
		const std::vector<BackBoneSite> &chain) {
	std::string positionmask;
	if (EnergyControls::getparameterset(parasetname_).keydefined("PositionMask"))
		EnergyControls::getnamedsetpar(parasetname_,"PositionMask", &positionmask);
	bool loopmask = (positionmask == "MaskLoopPositions");
	if (loopmask) {
		std::cout << "Loops will be  excluded from Energy calculations"
				<< std::endl;
	}
	positionmask_.clear();
	positionmask_.resize(chain.size(), false);
	int posi = 0;
	for (auto &s : chain) {
		if (loopmask) {
			char sscode = s.sscode;
			if (!(sscode == 'H' || sscode == 'E' || sscode == 'm'
					|| sscode == 'd'))
				positionmask_[posi] = true;
		}
		positionmask_[posi] = s.isgap || positionmask_[posi];
		++posi;
	}
}
ChainEnergyControl NSPproteinrep::prepareenergycontrol(
		std::vector<BackBoneSite> *chain,const std::string & parasetname) {
	ChainEnergyControl ce;
	int sspb;
	assert(EnergyControls::getparameterset(parasetname).initialized);
	if (EnergyControls::getparameterset(parasetname).keydefined("SSasPBType"))
		EnergyControls::getnamedsetpar(parasetname,"SSasPBType", &sspb);
	if (sspb != 0)
		ce.ssaspbtype() = true;
	int ignoreresname;
	if (EnergyControls::getparameterset(parasetname).keydefined("IgnoreResName"))
			EnergyControls::getnamedsetpar(parasetname,"IgnoreResName", &ignoreresname);
	if(ignoreresname != 0) ce.phipsi_ignoreresname()=true;
	std::string & origsscode = ce.origsscode();
	for (auto &s : *chain) {
		origsscode.push_back(s.sscode);
	}
	ce.setpositionmask(*chain);
	NSPpdbstatistics::ProteinBlock::setpbtypes(*chain, ce.ssaspbtype());
	std::string pbseq;
	for (auto &s : *chain) {
		pbseq.push_back(s.sscode);
	}
	std::string &refpbseq = ce.refpbseq();
	if (EnergyControls::getparameterset(parasetname).keydefined("RefPBSeq"))
		EnergyControls::getnamedsetpar(parasetname,"RefPBSeq", &refpbseq);
	if (refpbseq == "FromConf") {
		refpbseq.clear();
		refpbseq.resize(chain->size(), 'x');
		for (int i = 0; i < chain->size(); ++i) {
			char pbtype = chain->at(i).sscode;
			if (pbtype == 'm' || pbtype == 'd' || pbtype == 'H'
					|| pbtype == 'E') {
				refpbseq[i] = pbtype;
			}
		}
		std::cout << "starting sscode:  " << origsscode << std::endl;
		std::cout << "pbseq  from crd:  " << pbseq << std::endl;
		std::cout << "refpbseq for ene: " << refpbseq << std::endl;
	}

	for (int i = 0; i < chain->size(); ++i) {
		BackBoneSite &s = chain->at(i);
		s.resseq = i;
		if (i < chain->size() - 1) {
			if (chain->at(i + 1).resname == "PRO"
					|| chain->at(i + 1).resname == "CISPRO") {
				if (s.resname != "GLY" && s.resname != "PRO"
						&& s.resname != "CISPRO")
					s.resname = "PPR";
			}
		}
		if (i > 0) {
			if (chain->at(i - 1).nextpepcis())
				s.resname = "CISPRO";
		}
	}
//	EnergyControls::getnamedsetpar(parasetname,"PackingSeparation",&ce.blockpackingsep_);
	std::vector<double> & weights = ce.eweights();
	EnergyControls::getnamedsetpar(parasetname,"PhiPsiWeight", &weights[EnergyComponents::PHIPSI]);
	EnergyControls::getnamedsetpar(parasetname,"ClashWeight", &weights[EnergyComponents::CLASH]);
	EnergyControls::getnamedsetpar(parasetname,"PBLocalWeight",
			&weights[EnergyComponents::BLOCKLOCAL]);
	EnergyControls::getnamedsetpar(parasetname,"PBPackingWeight",
			&weights[EnergyComponents::BLOCKPACKING]);
	return ce;
}
double PhiPsiTerm::energy(double phi, double psi,
		const std::string &resname) const {
	const typename NSPpdbstatistics::PhiPsiDistr *distr;
	if(ignoreresname_)
		distr=&(NSPpdbstatistics::PhiPsiDistr::mixcoildistr());
	else if (resname == "PPR")
		distr = &(NSPpdbstatistics::PhiPsiDistr::phipsidistr("ALA", "PRO"));
	else
		distr = &(NSPpdbstatistics::PhiPsiDistr::phipsidistr(resname));
	return distr->statisticalenergy(phi, psi);
}
double PBlockTerm::energy(
		std::vector<BackBoneSite>::const_iterator iter) const {
	std::vector<double> torsions;
	int wh = windowwidth_ / 2;
	auto iter_begin = iter - wh;
	auto iter_end = iter_begin + windowwidth_;
	std::vector<double> pbtorsions;
	for (auto it = iter_begin; it != iter_end; ++it) {
		if (it == iter_begin) {
			if (it->omiga() > -90.0 && it->omiga() < 90.0)
				return 0.0; //cispep
			pbtorsions.push_back(it->psi());
		} else if (it == (iter_begin + windowwidth_ - 1))
			pbtorsions.push_back(it->phi());
		else {
			if (it->omiga() > -90.0 && it->omiga() < 90.0)
				return 0.0; //cispep
			pbtorsions.push_back(it->phi());
			pbtorsions.push_back(it->psi());
		}
	}
//	std::cout << " todo calculate energy for residue " << iter->resseq
//			<< std::endl;
	return energy(iter->resseq, iter->sscode, pbtorsions);

}
double PBlockTerm::energy(int resseq, char pbtype,
		const std::vector<double> &pbtorsions) const {
//	std::cout << " pbtype: " << pbtype << " torsions:";
//	for (auto t : pbtorsions)
//		std::cout << " " << t;
//	std::cout << std::endl;
	if (pbtype == 't')
		return 0.0;
	double wrongpbene = 0.01;
	if (!refpbseq_->empty()) {
		char expectedpb = (*refpbseq_)[resseq];
		if (expectedpb != 'x' && expectedpb != pbtype) {
			std::vector<double> dev =
					NSPpdbstatistics::ProteinBlock::deviations(expectedpb,
							pbtorsions);
			double dev2 = 0.0;
			for (auto d : dev)
				dev2 += d * d;
			dev2 /= 8.0;
			return wrongpbene * dev2;
		}
	}
	return NSPpdbstatistics::NN_PepScorer::energy(pbtype, pbtorsions);
}
double LocalBackBoneEnergy::totalenergy(const std::vector<BackBoneSite> & chain,
		ChainEnergyControl *ce, int ignorehead, int ignoretail) {
	if (terms_.empty())
		return 0.0;
	double ene = 0.0;
	auto iterend = chain.begin() + (chain.size() - ignoretail);
	std::vector<double> ec(terms_.size(), 0.0);
	for (auto iter = chain.begin() + ignorehead; iter != iterend; ++iter) {
		if (ce->positionmasked(iter->resseq))
			continue;
		int tidx = 0;
		for (auto &term : terms_) {
			int ww = term->windowwidth();
			int wh = ww / 2;
			if (iter - chain.begin() < wh || chain.end() - iter + wh < ww)
				continue;
			double e = term->energy(iter) * ce->eweight(comptypes_[tidx]);
			ec[tidx++] += e;
		}
	}
	int tidx = 0;
	for (auto e : ec) {
		ene += e;
		ce->energy(comptypes_[tidx++]) += e;
	}
	return ene;
}
double LocalBackBoneEnergy::partialE(const std::vector<BackBoneSite> & chain,
		int movestart, int moveend, const std::vector<BackBoneSite> &moved,
		ChainEnergyControl *ce) {
	if (terms_.empty())
		return 0.0;
	double ene = 0.0;
	int ww = wwmax_;
	if (movestart < moveend) {
		int beginposi = movestart - ww + 1;
		if (beginposi < 0)
			beginposi = 0;
		int endposi = moveend + ww - 1;
		if (endposi > chain.size())
			endposi = chain.size();
		std::vector<BackBoneSite> newpart(endposi - beginposi);
		if (movestart > beginposi)
			std::copy(chain.begin() + beginposi, chain.begin() + movestart,
					newpart.begin());
		std::copy(moved.begin(), moved.end(),
				newpart.begin() + movestart - beginposi);
		if (endposi > moveend)
			std::copy(chain.begin() + moveend, chain.begin() + endposi,
					newpart.begin() + moveend - beginposi);
		int ignorehead = wwh_ + (movestart - ww + 1 - beginposi);
		int ignoretail = wwh_ + (endposi - moveend - ww + 1);
		if (ignorehead < 0)
			ignorehead = 0;
		if (ignoretail < 0)
			ignoretail = 0;
		NSPpdbstatistics::ProteinBlock::setpbtypes(newpart, ce->ssaspbtype());
		ene = totalenergy(newpart, ce, ignorehead, ignoretail);
	} else {
		int beginposi = movestart - ww + 1;
		if (beginposi < 0)
			beginposi = 0;
		std::vector<BackBoneSite> newpart(chain.size() - beginposi);
		if (movestart > beginposi)
			std::copy(chain.begin() + beginposi, chain.begin() + movestart,
					newpart.begin());
		std::copy(moved.begin(), moved.begin() + chain.size() - movestart,
				newpart.begin() + movestart - beginposi);
		int ignorehead = wwh_ + (movestart - ww + 1 - beginposi);
		if (ignorehead < 0)
			ignorehead = 0;
		NSPpdbstatistics::ProteinBlock::setpbtypes(newpart, ce->ssaspbtype());
		ene = totalenergy(newpart, ce, ignorehead, 0);
		if (moveend > 0) {
			int endposi = moveend + ww - 1;
			if (endposi > chain.size())
				endposi = chain.size();
			std::vector<BackBoneSite> newpart(endposi);
			std::copy(moved.begin() + chain.size() - movestart, moved.end(),
					newpart.begin());
			if (endposi > moveend)
				std::copy(chain.begin() + moveend, chain.begin() + endposi,
						newpart.begin() + moveend);
			int ignoretail = wwh_ + (endposi - moveend - ww + 1);
			if (ignoretail < 0)
				ignoretail = 0;
			NSPpdbstatistics::ProteinBlock::setpbtypes(newpart,
					ce->ssaspbtype());
			ene += totalenergy(newpart, ce, 0, ignoretail);
		}
	}
	return ene;
}
double LocalBackBoneEnergy::partialE(const std::vector<BackBoneSite> & chain,
		int movestart, int moveend, ChainEnergyControl *ce) {
	if (terms_.empty())
		return 0.0;
	double ene = 0.0;
	std::vector<double> ec(terms_.size(), 0.0);
	int ww = wwmax_;
	if (moveend < movestart) {
		int beginposi = movestart - ww + 1;
		if (beginposi < 0)
			beginposi = 0;
		int ignorehead = wwh_ + (movestart - ww + 1 - beginposi);
		if (ignorehead < 0)
			ignorehead = 0;
		for (auto iter = chain.begin() + beginposi + ignorehead;
				iter != chain.end(); ++iter) {
			if (ce->positionmasked(iter->resseq))
				continue;
			int tidx = 0;
			for (auto &term : terms_) {
				int ww = term->windowwidth();
				int wh = ww / 2;
				if (iter - chain.begin() - beginposi < wh
						|| chain.end() - iter + wh < ww)
					continue;
				double e = term->energy(iter) * ce->eweight(comptypes_[tidx]);
				ec[tidx++] += e;
			}
		}
		if (moveend > 0) {
			int endposi = moveend + ww - 1;
			if (endposi > chain.size())
				endposi = chain.size();
			int ignoretail = wwh_ + (endposi - moveend - ww + 1);
			if (ignoretail < 0)
				ignoretail = 0;
			for (auto iter = chain.begin();
					iter != chain.begin() + endposi - ignoretail; ++iter) {
				if (ce->positionmasked(iter->resseq))
					continue;
				int tidx = 0;
				for (auto &term : terms_) {
					int ww = term->windowwidth();
					int wh = ww / 2;
					if (iter - chain.begin() < wh
							|| chain.begin() + endposi - iter + wh < ww)
						continue;
					double e = term->energy(iter)
							* ce->eweight(comptypes_[tidx]);
					ec[tidx++] += e;
				}
			}
		}
	} else {
		int beginposi = movestart - ww + 1;
		if (beginposi < 0)
			beginposi = 0;
		int endposi = moveend + ww - 1;
		if (endposi > chain.size())
			endposi = chain.size();
		int ignorehead = wwh_ + (movestart - ww + 1 - beginposi);
		if (ignorehead < 0)
			ignorehead = 0;
		int ignoretail = wwh_ + (endposi - moveend - ww + 1);
		if (ignoretail < 0)
			ignoretail = 0;
		for (auto iter = chain.begin() + beginposi + ignorehead;
				iter != chain.begin() + endposi - ignoretail; ++iter) {
//			if(iter->resseq <0 || iter->resseq >= chain.size()) {
//				std::cout <<beginposi << " "<<endposi<<" " << ignorehead<<" " << ignoretail <<std::endl;
//				abort();
//			}
			if (ce->positionmasked(iter->resseq))
				continue;
			int tidx = 0;
			for (auto &term : terms_) {
				int ww = term->windowwidth();
				int wh = ww / 2;
				if (iter - chain.begin() - beginposi < wh
						|| chain.begin() + endposi - iter + wh < ww)
					continue;
				double e = term->energy(iter) * ce->eweight(comptypes_[tidx]);
				ec[tidx++] += e;
			}
		}
	}
	int tidx = 0;
	for (auto e : ec) {
		ene += e;
		ce->energy(comptypes_[tidx++]) += e;
	}
	return ene;
}
/*static void intersitedistmatrix(const BackBoneSite &s1, const BackBoneSite &s2,
 std::vector<std::vector<double>> *dist2matrix) {
 dist2matrix->clear();
 std::vector<NSPgeometry::XYZ> crd1, crd2;
 s1.getcrd(crd1);
 s2.getcrd(crd2);
 for (int i = 0; i < crd1.size(); ++i) {
 dist2matrix->push_back(std::vector<double>());
 std::vector<double> &row = dist2matrix->back();
 for (int j = 0; j < crd2.size(); ++j)
 row.push_back((crd1[i] - crd2[j]).squarednorm());
 }
 }*/
double StericClashTerm::energy(std::vector<BackBoneSite>::const_iterator iter1,
		std::vector<BackBoneSite>::const_iterator iter2,
		double *cadist2) const {
	int sep = iter1->resseq - iter2->resseq;
	if (sep < minsep_ && sep > -minsep_)
		return 0.0;
	if (*cadist2 <= 0)
		*cadist2 = (iter1->cacrd() - iter2->cacrd()).squarednorm();
	if (*cadist2 > cacutoff2_)
		return 0.0;
//	if(print) {
//		std::cout << "determine clash between " << iter1->resseq << " v.s. "
//			<< iter2->resseq << std::endl;
//	}
//	return 1.0;
	bool clashed = false;
	for (int i = 0; i < 4; ++i) {
		NSPgeometry::XYZ a1 = iter1->getcrd(3 * i + BackBoneSite::NCRD);
		for (int j = 0; j < 4; ++j) {
			NSPgeometry::XYZ a2 = iter2->getcrd(3 * j + BackBoneSite::NCRD);
			double d2 = (a1 - a2).squarednorm();
			if (d2 <= 8.41) {
				if ((i == 0 && j == 3) || (i == 3 && j == 0)) {
					if (d2 <= 5.76)
						clashed = true;
				} else
					clashed = true;
			}
			if (clashed) {
//			std::cout <<"Pairs clashed: "<<i <<" "<<j <<" "<<d2<<std::endl;
//				std::cout<< iter1->toString();
//				std::cout <<iter2->toString();
				return clashenergy_;
			}
		}
	}
	return 0.0;
}
double StericClashTerm::energy(const BackBoneSite &s1, const BackBoneSite &s2,
		double *cadist2) const {
	if (*cadist2 <= 0)
		*cadist2 = (s1.cacrd() - s2.cacrd()).squarednorm();
	if (*cadist2 > cacutoff2_)
		return 0.0;
	bool clashed = false;
	for (int i = 0; i < 4; ++i) {
		NSPgeometry::XYZ a1 = s1.getcrd(3 * i + BackBoneSite::NCRD);
		for (int j = 0; j < 4; ++j) {
			NSPgeometry::XYZ a2 = s2.getcrd(3 * j + BackBoneSite::NCRD);
			double d2 = (a1 - a2).squarednorm();
			if (d2 <= 8.41) {
				if ((i == 0 && j == 3) || (i == 3 && j == 0)) {
					if (d2 <= 5.76)
						clashed = true;
				} else
					clashed = true;
			}
			if (clashed) {
				return clashenergy_;
			}
		}
	}
	return 0.0;
}
double BlockPackingTerm::energy(std::vector<BackBoneSite>::const_iterator iter1,
		std::vector<BackBoneSite>::const_iterator iter2,
		double *cadist2) const {
	int sep = iter1->resseq - iter2->resseq;
	if (sep < minsep_ && sep > -minsep_)
		return 0.0;
	if (*cadist2 <= 0)
		*cadist2 = (iter1->cacrd() - iter2->cacrd()).squarednorm();
	if (*cadist2 > cacutoff2_)
		return 0.0;
	if (iter1->sscode == 't' || iter2->sscode == 't')
		return 0.0;
//todo real energy
//	std::cout << "determine blockpacking between " << iter1->resseq << "-"
//			<< iter1->sscode << " v.s. " << iter2->resseq << " "
//			<< iter2->sscode << std::endl;
	double e;
	if (iter1->resseq < iter2->resseq)
		e = NSPpdbstatistics::tetrapairenergy(*iter1, *iter2);
	else
		e = NSPpdbstatistics::tetrapairenergy(*iter2, *iter1);
	/*	if(e>10.0) {
	 std::cout<< iter1->resid<<iter1->resname<<"-"<<iter2->resid<<iter2->resname<<": ";
	 NSPpdbstatistics::TetraGeom geom(*iter1,*iter2);
	 std::cout <<geom.pbtypes().first<<" "<< geom.pbtypes().second
	 <<" "<<geom.orient().first<<" "<<geom.orient().second
	 <<" "<<geom.rmin() <<std::endl;
	 }*/
	return e;
}
double BlockPackingTerm::energy(const BackBoneSite &s1, const BackBoneSite &s2,
		double *cadist2) const {
	if (*cadist2 <= 0)
		*cadist2 = (s1.cacrd() - s2.cacrd()).squarednorm();
	if (*cadist2 > cacutoff2_)
		return 0.0;
	if (s1.sscode == 't' || s2.sscode == 't')
		return 0.0;
	double e;
	if (s1.resseq < s2.resseq)
		e = NSPpdbstatistics::tetrapairenergy(s1, s2);
	else
		e = NSPpdbstatistics::tetrapairenergy(s2, s1);
	return e;
}
double BackBonePackingEnergy::totalenergy(
		const std::vector<BackBoneSite> & chain, ChainEnergyControl *ce,
		int ignorehead, int ignoretail) {
	if (terms_.empty())
		return 0.0;
	double ene = 0.0;
	std::vector<double> ec(terms_.size(), 0.0);
	auto iterend = chain.begin() + (chain.size() - ignoretail);
	for (auto iter1 = chain.begin() + ignorehead; iter1 != iterend; ++iter1) {
		if (ce->positionmasked(iter1->resseq))
			continue;
		for (auto iter2 = iter1 + 1; iter2 != iterend; ++iter2) {
			if (ce->positionmasked(iter2->resseq))
				continue;
			int tidx = 0;
			double cadist2 = -1.0;
			for (auto &term : terms_) {
				int ww = term->windowwidth();
				int wh = ww / 2;
				if (iter1 - chain.begin() < wh || chain.end() - iter1 + wh < ww)
					continue;
				if (iter2 - chain.begin() < wh || chain.end() - iter2 + wh < ww)
					continue;
				double e = term->energy(iter1, iter2, &cadist2)
						* ce->eweight(comptypes_[tidx]);
				ec[tidx++] += e;
			}
		}
	}
	int tidx = 0;
	for (auto e : ec) {
		ene += e;
		ce->energy(comptypes_[tidx++]) += e;
	}
//	std::cout << "Packing total ENE: " << ene << std::endl;
	return ene;
}
double BackBonePackingEnergy::partialE(const std::vector<BackBoneSite> & chain,
		int movestart, int moveend, const std::vector<BackBoneSite> &moved,
		ChainEnergyControl *ce) {
	if (terms_.empty())
		return 0.0;
	double ene = 0.0;
	int ww = wwmax_;
	std::vector<double> ec(terms_.size(), 0.0);
	if (movestart < moveend) {
		int beginposi = movestart - ww + 1;
		if (beginposi < 0)
			beginposi = 0;
		int endposi = moveend + ww - 1;
		if (endposi > chain.size())
			endposi = chain.size();
		std::vector<BackBoneSite> newpart(endposi - beginposi);
		if (movestart > beginposi)
			std::copy(chain.begin() + beginposi, chain.begin() + movestart,
					newpart.begin());
		std::copy(moved.begin(), moved.end(),
				newpart.begin() + movestart - beginposi);
		if (endposi > moveend)
			std::copy(chain.begin() + moveend, chain.begin() + endposi,
					newpart.begin() + moveend - beginposi);
		int ignorehead = wwh_ + (movestart - ww + 1 - beginposi);
		if (ignorehead < 0)
			ignorehead = 0;
		int ignoretail = wwh_ + (endposi - moveend - ww + 1);
		if (ignoretail < 0)
			ignoretail = 0;
		NSPpdbstatistics::ProteinBlock::setpbtypes(newpart, ce->ssaspbtype());
		ene += totalenergy(newpart, ce, ignorehead, ignoretail);
		for (auto iter1 = chain.begin(); iter1 != chain.end(); ++iter1) {
			if (ce->positionmasked(iter1->resseq))
				continue;
			if (iter1 - chain.begin() >= beginposi + ignorehead
					&& iter1 - chain.begin() < endposi - ignoretail)
				continue;
			for (auto iter2 = newpart.begin() + ignorehead;
					iter2 != newpart.end() - ignoretail; ++iter2) {
				if (ce->positionmasked(iter2->resseq))
					continue;
				int tidx = 0;
				double cadist2 = -1.0;
				for (auto &term : terms_) {
					int ww = term->windowwidth();
					int wh = ww / 2;
					if (iter1 - chain.begin() < wh
							|| chain.end() - iter1 + wh < ww)
						continue;
					if (iter2 - newpart.begin() < wh
							|| newpart.end() - iter2 + wh < ww)
						continue;
					double e = term->energy(iter1, iter2, &cadist2)
							* ce->eweight(comptypes_[tidx]);
					ec[tidx++] += e;
				}
			}
		}
	} else {
		int beginposi = movestart - ww + 1;
		if (beginposi < 0)
			beginposi = 0;
		int endposi = moveend + ww - 1;
		if (endposi > chain.size())
			endposi = chain.size();
		std::vector<BackBoneSite> newpart(chain.size() - beginposi);
		if (movestart > beginposi)
			std::copy(chain.begin() + beginposi, chain.begin() + movestart,
					newpart.begin());
		std::copy(moved.begin(), moved.begin() + chain.size() - movestart,
				newpart.begin() + movestart - beginposi);
		int ignorehead = wwh_ + (movestart - ww + 1 - beginposi);
		if (ignorehead < 0)
			ignorehead = 0;
		int ignoretail = wwh_ + (endposi - moveend - ww + 1);
		if (ignoretail < 0)
			ignoretail = 0;
		NSPpdbstatistics::ProteinBlock::setpbtypes(newpart, ce->ssaspbtype());
		ene += totalenergy(newpart, ce, ignorehead, 0);
		auto it1start = chain.begin() + endposi - ignoretail;
		if (moveend == 0)
			it1start = chain.begin();
		for (auto iter1 = it1start; iter1 != chain.end(); ++iter1) {
			if (ce->positionmasked(iter1->resseq))
				continue;
			if (iter1 - chain.begin() >= beginposi + ignorehead)
				continue;
			for (auto iter2 = newpart.begin() + ignorehead;
					iter2 != newpart.end(); ++iter2) {
				if (ce->positionmasked(iter2->resseq))
					continue;
				int tidx = 0;
				double cadist2 = -1.0;
				for (auto &term : terms_) {
					int ww = term->windowwidth();
					int wh = ww / 2;
					if (iter1 - chain.begin() < wh
							|| chain.end() - iter1 + wh < ww)
						continue;
					if (iter2 - newpart.begin() < wh
							|| newpart.end() - iter2 + wh < ww)
						continue;
					double e = term->energy(iter1, iter2, &cadist2)
							* ce->eweight(comptypes_[tidx]);
					ec[tidx++] += e;
				}
			}
		}
		if (moveend > 0) {
			std::vector<BackBoneSite> newpart2(endposi);
			std::copy(moved.begin() + chain.size() - movestart, moved.end(),
					newpart2.begin());
			if (endposi > moveend)
				std::copy(chain.begin() + moveend, chain.begin() + endposi,
						newpart2.begin() + moveend);
			NSPpdbstatistics::ProteinBlock::setpbtypes(newpart2,
					ce->ssaspbtype());
			ene += totalenergy(newpart2, ce, 0, ignoretail);
			for (auto iter1 = chain.begin();
					iter1 != chain.begin() + beginposi + ignorehead; ++iter1) {
				if (ce->positionmasked(iter1->resseq))
					continue;
				if (iter1 - chain.begin() < endposi - ignoretail)
					continue;
				for (auto iter2 = newpart2.begin();
						iter2 != newpart2.end() - ignoretail; ++iter2) {
					if (ce->positionmasked(iter2->resseq))
						continue;
					int tidx = 0;
					double cadist2 = -1.0;
					for (auto &term : terms_) {
						int ww = term->windowwidth();
						int wh = ww / 2;
						if (iter1 - chain.begin() < wh
								|| chain.end() - iter1 + wh < ww)
							continue;
						if (iter2 - newpart2.begin() < wh
								|| newpart2.end() - iter2 + wh < ww)
							continue;
						double e = term->energy(iter1, iter2, &cadist2)
								* ce->eweight(comptypes_[tidx]);
						ec[tidx++] += e;
					}
				}
			}
			for (auto iter1 = newpart.begin() + ignorehead;
					iter1 != newpart.end(); ++iter1) {
				if (ce->positionmasked(iter1->resseq))
					continue;
				for (auto iter2 = newpart2.begin();
						iter2 != newpart2.end() - ignoretail; ++iter2) {
					if (ce->positionmasked(iter2->resseq))
						continue;
					int tidx = 0;
					double cadist2 = -1.0;
					for (auto &term : terms_) {
						int ww = term->windowwidth();
						int wh = ww / 2;
						if (iter1 - newpart.begin() < wh
								|| newpart.end() - iter1 + wh < ww)
							continue;
						if (iter2 - newpart2.begin() < wh
								|| newpart2.end() - iter2 + wh < ww)
							continue;
						double e = term->energy(iter1, iter2, &cadist2)
								* ce->eweight(comptypes_[tidx]);
						ec[tidx++] += e;
					}
				}
			}
		}
	}
//	std::cout << "Packing partial ENE: " << ene << std::endl;
	int tidx = 0;
	for (auto e : ec) {
		ene += e;
		ce->energy(comptypes_[tidx++]) += e;
	}
	return ene;
}
double BackBonePackingEnergy::partialE(const std::vector<BackBoneSite> & chain,
		int movestart, int moveend, ChainEnergyControl *ce) {
	if (terms_.empty())
		return 0.0;
	double ene = 0.0;
	int ww = wwmax_;
	int wh = wwh_;
	std::vector<double> ec(terms_.size(), 0.0);
	for (auto iter1 = chain.begin(); iter1 != chain.end(); ++iter1) {
		if (ce->positionmasked(iter1->resseq))
			continue;
		int beginposi = movestart - ww + 1;
		if (beginposi < 0)
			beginposi = 0;
		int endposi = moveend + ww - 1;
		if (endposi > chain.size())
			endposi = chain.size();
		int ignorehead = wwh_ + (movestart - beginposi - ww + 1);
		if (ignorehead < 0)
			ignorehead = 0;
		int ignoretail = wwh_ + (endposi - moveend - ww + 1);
		if (ignoretail < 0)
			ignoretail = 0;
		if (moveend < movestart) {
			bool changed1 = (iter1 - chain.begin() >= beginposi + ignorehead)
					|| (iter1 - chain.begin() < endposi - ignoretail);
			for (auto iter2 = chain.begin() + beginposi + ignorehead;
					iter2 != chain.end(); ++iter2) {
				if (changed1 && iter1 - iter2 >= 0)
					continue;
				if (ce->positionmasked(iter2->resseq))
					continue;
				int tidx = 0;
				double cadist2 = -1.0;
				for (auto &term : terms_) {
					int ww = term->windowwidth();
					int wh = ww / 2;
					if (iter1 - chain.begin() < wh
							|| chain.end() - iter1 + wh < ww)
						continue;
					if (iter2 - chain.begin() - beginposi < wh
							|| chain.end() - iter2 + wh < ww)
						continue;
					double e = term->energy(iter1, iter2, &cadist2)
							* ce->eweight(comptypes_[tidx]);
					ec[tidx++] += e;
				}
			}
			if (moveend > 0) {
				for (auto iter2 = chain.begin();
						iter2 != chain.begin() + endposi - ignoretail;
						++iter2) {
					if (changed1 && iter1 - iter2 >= 0)
						continue;
					if (ce->positionmasked(iter2->resseq))
						continue;
					int tidx = 0;
					double cadist2 = -1.0;
					for (auto &term : terms_) {
						int ww = term->windowwidth();
						int wh = ww / 2;
						if (iter1 - chain.begin() < wh
								|| chain.end() - iter1 + wh < ww)
							continue;
						if (iter2 - chain.begin() < wh
								|| chain.begin() + endposi - iter2 + wh < ww)
							continue;
						double e = term->energy(iter1, iter2, &cadist2)
								* ce->eweight(comptypes_[tidx]);
						ec[tidx++] += e;
					}
				}
			}
		} else {
			bool changed1 = iter1 - chain.begin() >= beginposi + ignorehead
					&& iter1 - chain.begin() < endposi - ignoretail;
			for (auto iter2 = chain.begin() + beginposi + ignorehead;
					iter2 != chain.begin() + endposi - ignoretail; ++iter2) {
				if (changed1 && iter1 - iter2 >= 0)
					continue;
				if (ce->positionmasked(iter2->resseq))
					continue;
				int tidx = 0;
				double cadist2 = -1.0;
				for (auto &term : terms_) {
					int ww = term->windowwidth();
					int wh = ww / 2;
					if (iter1 - chain.begin() < wh
							|| chain.end() - iter1 + wh < ww)
						continue;
					if (iter2 - chain.begin() - beginposi < wh
							|| chain.begin() + endposi - iter2 + wh < ww)
						continue;
					double e = term->energy(iter1, iter2, &cadist2)
							* ce->eweight(comptypes_[tidx]);
					ec[tidx++] += e;
				}
			}
		}
	}
	int tidx = 0;
	for (auto e : ec) {
		ene += e;
		ce->energy(comptypes_[tidx++]) += e;
	}
	return ene;
}
double BackBoneEnergy::partialE(const std::vector<BackBoneSite> &chain, int partstart,
		int partend,ChainEnergyControl *ce) {
	lene_.setrefpbseq(&(ce->refpbseq()));
//	if(moveend==0) print=true;
	EnergyComponents esave = ce->ecomp();

	EnergyComponents ec0;
	ce->ecomp() = ec0;
	double ene0 = lene_.partialE(chain, partstart, partend, ce);
	ene0 += pene_.partialE(chain, partstart, partend, ce);
	ce->ecomp().energies[EnergyComponents::TOTAL]=ene0;
	return ene0;
}
double BackBoneEnergy::deltaE(BackBoneMoves &moves, int mvidx,
		ChainEnergyControl *ce, EnergyComponents *ecomp0,
		EnergyComponents *ecomp1) {
	int movestart = moves.startposi;
	int moveend = moves.endposi;
	lene_.setrefpbseq(&(ce->refpbseq()));
//	if(moveend==0) print=true;
	EnergyComponents esave = ce->ecomp();
	if (!moves.energy0set) {
//		std::cout << "ref calc" << std::endl;
		EnergyComponents ec0;
		ce->ecomp() = ec0;
		double ene0 = lene_.partialE(*moves.hostchain, movestart, moveend, ce);
		ene0 += pene_.partialE(*moves.hostchain, movestart, moveend, ce);
		moves.energy0 = ene0;
		*ecomp0 = ce->ecomp();
		ecomp0->energies[EnergyComponents::TOTAL] = ene0;
		moves.energy0set = true;
//		std::cout <<"refNterms: " <<ecomp0->energies[EnergyComponents::CLASH]<<std::endl;
	}
//	std::cout << "moved_calc" << std::endl;
	EnergyComponents newec;
	ce->ecomp() = newec;
	double enenew = lene_.partialE(*moves.hostchain, movestart, moveend,
			*moves.movedloop(mvidx), ce);
	enenew += pene_.partialE(*moves.hostchain, movestart, moveend,
			*moves.movedloop(mvidx), ce);
	*ecomp1 = ce->ecomp();
	ecomp1->energies[EnergyComponents::TOTAL] = enenew;
//	std::cout <<"movedNterms: " <<ecomp1->energies[EnergyComponents::CLASH]<<std::endl;
//	if(ecomp1->energies[EnergyComponents::CLASH] != ecomp0->energies[EnergyComponents::CLASH]){
//		std::cout <<moves.startposi <<" " <<moves.endposi<<std::endl;
//		abort();
//	}
	ce->ecomp() = esave;
	return enenew - moves.energy0;
}
static bool negligible(double w) {
	static const double EPSINON = 1.e-16;
	if (w < EPSINON && w > -EPSINON)
		return true;
	return false;
}
LocalBackBoneEnergy NSPproteinrep::makelocalbackboneenergy(
		const std::vector<double> &weights,bool phipsi_ignoreresname) {
	LocalBackBoneEnergy res;
	std::shared_ptr<BackBoneLocalTerm> term1 = std::shared_ptr
			< BackBoneLocalTerm > (new PhiPsiTerm(phipsi_ignoreresname));
	if (!negligible(weights[EnergyComponents::PHIPSI])){
		res.addterm(term1, EnergyComponents::PHIPSI);
	}
	std::shared_ptr<BackBoneLocalTerm> term2 = std::shared_ptr
			< BackBoneLocalTerm > (new PBlockTerm());
	if (!negligible(weights[EnergyComponents::BLOCKLOCAL]))
		res.addterm(term2, EnergyComponents::BLOCKLOCAL);
	return res;
}

BackBonePackingEnergy NSPproteinrep::makebackbonepackingenergy(
		const std::vector<double> &weights) {
	BackBonePackingEnergy res;
	std::shared_ptr<BackBonePackingTerm> term1 = std::shared_ptr
			< BackBonePackingTerm > (new StericClashTerm());
	if (!negligible(weights[EnergyComponents::CLASH]))
		res.addterm(term1, EnergyComponents::CLASH);
	auto & eparset=EnergyControls::getparameterset();
	if(!eparset.initialized) initenergycontrols();
	int packingsep;
	eparset.getval("PackingSeparation",&packingsep);
	std::shared_ptr<BackBonePackingTerm> term2 = std::shared_ptr
			< BackBonePackingTerm > (new BlockPackingTerm(packingsep));
	if (!negligible(weights[EnergyComponents::BLOCKPACKING]))
		res.addterm(term2, EnergyComponents::BLOCKPACKING);
	return res;
}
double BackBonePackingEnergy::totalenergy(const ChainPack &cpck,
		EnergyComponents *ecomp) {
	int idx = 0;
	double etotal = 0.0;
	for (auto &chain : cpck) {
		double e = totalenergy(chain, &(cpck.econtrol(idx)));
		cpck.energy(idx, idx) += e;
		etotal += e;
		for (int tidx = 0; tidx < terms_.size(); ++tidx) {
			ecomp->energies[comptypes_[tidx]] += cpck.econtrol(idx).energy(
					comptypes_[tidx]);
		}
		++idx;
	}
	std::vector<double> &eweights = cpck.econtrol(0).eweights();
	std::vector<double> ec(terms_.size(), 0.0);
	for (int i = 0; i < cpck.size() - 1; ++i) {
		for (auto &s1 : cpck.at(i)) {
			for (int j = i + 1; j < cpck.size(); ++j) {
				double em = 0.0;
				for (auto &s2 : cpck.at(j)) {
					int tidx = 0;
					double cadist2 = -1.0;
					for (auto &term : terms_) {
						double e = term->energy(s1, s2, &cadist2)
								* eweights[comptypes_[tidx]];
						em += e;
						ec[tidx++] += e;
					}
				}
				cpck.energy(i, j) += em;
			}
		}
	}
	int tidx = 0;
	for (auto e : ec) {
		etotal += e;
		ecomp->energies[comptypes_[tidx++]] += e;
	}
	return etotal;
}
void BackBonePackingEnergy::partialE(const ChainPack &cpck,
		ChainPackMoves *moves, int mvidx) {
	for (auto c : moves->moveset) {
		if (!moves->localmove[c])
			continue;
		double e = totalenergy(moves->movedchain(c, mvidx), &cpck.econtrol(c));
		moves->energy(c, c, mvidx) += e;
	}
	std::vector<double> &eweights = cpck.econtrol(0).eweights();
	for (int i = 0; i < cpck.size(); ++i) {
		if (moves->moveset.find(i) != moves->moveset.end())
			continue;
		for (auto &s1 : cpck.at(i)) {
			for (auto c : moves->moveset) {
				double em = 0.0;
				for (BackBoneSite &s2 : moves->movedchain(c, mvidx)) {
					int tidx = 0;
					double cadist2 = -1.0;
					for (auto &term : terms_) {
						if (c > i) {
							double e = term->energy(s1, s2, &cadist2)
									* eweights[comptypes_[tidx++]];
							em += e;
						} else {
							double e = term->energy(s2, s1, &cadist2)
									* eweights[comptypes_[tidx++]];
							em += e;
						}
					}
				}
				moves->energy(i, c, mvidx) += em;
			}
		}
	}
}
double LocalBackBoneEnergy::totalenergy(const ChainPack & cpck,
		EnergyComponents *ecomp) {
	int idx = 0;
	double etotal = 0.0;
	for (auto &chain : cpck) {
		setrefpbseq(&(cpck.econtrol(idx).refpbseq()));
		double e = totalenergy(chain, &(cpck.econtrol(idx)));
		cpck.energy(idx, idx) += e;
		etotal += e;
		for (int tidx = 0; tidx < terms_.size(); ++tidx) {
			ecomp->energies[comptypes_[tidx]] += cpck.econtrol(idx).energy(
					comptypes_[tidx]);
		}
		++idx;
	}
	return etotal;
}
void LocalBackBoneEnergy::partialE(const ChainPack &cpck, ChainPackMoves *moves,
		int mvidx) {
	for (auto c : moves->moveset) {
		if (!moves->localmove[c])
			continue;
		setrefpbseq(&(cpck.econtrol(c).refpbseq()));
		double e = totalenergy(moves->movedchain(c, mvidx), &cpck.econtrol(c));
		moves->energy(c, c, mvidx) += e;
	}
}
double BackBoneEnergy::totalenergy(const ChainPack & cpck,
		EnergyComponents *ecomp) {
	ecomp->reset();
	for (int i = 0; i < cpck.size(); ++i) {
		cpck.econtrol(i).resetecomp();
		cpck.initematrix();
	}
	double etot = lene_.totalenergy(cpck, ecomp)
			+ pene_.totalenergy(cpck, ecomp);
	ecomp->energies[EnergyComponents::TOTAL] = etot;
	return etot;
}
double BackBoneEnergy::deltaE(const ChainPack & cpck, ChainPackMoves *moves,
		int mvidx) {
	moves->ematrices[mvidx].init(cpck.size());
	lene_.partialE(cpck, moves);
	pene_.partialE(cpck, moves);
	double de = 0.0;
	for (auto c : moves->moveset) {
		if (moves->localmove[c])
			de += moves->energy(c, c, mvidx) - cpck.energy(c, c);
		for (int j = 0; j < cpck.size(); ++j) {
			if (moves->moveset.find(j) != moves->moveset.end())
				continue;
			de += moves->energy(c, j, mvidx) - cpck.energy(c, j);
		}
	}
	return de;
}
