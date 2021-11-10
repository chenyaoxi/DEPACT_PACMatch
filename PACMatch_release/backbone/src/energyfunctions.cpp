/*
 * energyfunctions.cpp
 *
 *  Created on: 2017年4月25日
 *      Author: hyliu
 */
#include "backbone/mainchain.h"
#include "backbone/energyfunctions.h"
#include "backbone/rminsef.h"
#include "pdbstatistics/phipsidistr.h"
#include "backbone/torsionvectorscorer.h"
#include "dataio/datapaths.h"
using namespace NSPproteinrep;
using namespace NSPgeometry;
int LocalEne::fragmentlength { 5 };
double SiteAtomContact::energy(const BackBoneSite & bs1,
		const BackBoneSite &bs2, int sep) {
	if (sep < minsep_)
		return 0.0;
	NSPgeometry::XYZ r1 = bs1.getcrd(BackBoneSite::NCRD + atom1_ * 3);
	NSPgeometry::XYZ r2 = bs2.getcrd(BackBoneSite::NCRD + atom2_ * 3);
	double d2 = NSPgeometry::distance2(r1, r2);
	if (d2 < rcut2_)
		return 1.0;
	else
		return 0.0;
}

TorsionEne::TorsionEne() {
// read phipsi  distributions from data files
	std::string datapath = NSPdataio::datapath();
	NSPpdbstatistics::PhiPsiDistr::coildistr(datapath + "coilphipsi.dat");
	NSPpdbstatistics::PhiPsiDistr::glydistr(datapath + "glyphipsi.dat");
	NSPpdbstatistics::PhiPsiDistr::preprodistr(datapath + "preprophipsi.dat");
	NSPpdbstatistics::PhiPsiDistr::transprodistr(
			datapath + "transprophipsi.dat");
	NSPpdbstatistics::PhiPsiDistr::cisprodistr(datapath + "cisprophipsi.dat");
}

double TorsionEne::energy(const std::vector<const BackBoneSite *> &sites) {
	assert(sites.size() == LocalEne::fragmentlength);
	if (!defined(sites))
		return 0.0;
	int posim = LocalEne::fragmentlength / 2;
	const NSPpdbstatistics::PhiPsiDistr *distr;
	distr = &(NSPpdbstatistics::PhiPsiDistr::coildistr());
	std::string resname = sites[posim]->resname;
	if (resname == "GLY") {
		distr = &(NSPpdbstatistics::PhiPsiDistr::glydistr());
	} else if (resname == "PRO") {
		double omiga = 180.0;
		if (sites[posim - 1])
			omiga = sites[posim - 1]->omiga();
		if (omiga > -90.0 && omiga < 90.0)
			distr = &(NSPpdbstatistics::PhiPsiDistr::cisprodistr());
		else
			distr = &(NSPpdbstatistics::PhiPsiDistr::transprodistr());
	} else if (sites[posim + 1]) {
		if (sites[posim + 1]->resname == "PRO")
			distr = &(NSPpdbstatistics::PhiPsiDistr::preprodistr());
	}
	return distr->statisticalenergy(sites[posim]->phi(), sites[posim]->psi());
}

TorsionVecEne::TorsionVecEne() {
// read phipsi  distributions from data files
	std::string datapath = NSPdataio::datapath();
	std::vector<BackBoneSite> tmpsites;
	readbackbonesites(datapath + "tmplatesites.dat", tmpsites);
	TorsionVectorScorer::getinstance(&tmpsites);
}
double TorsionVecEne::energy(const std::vector<const BackBoneSite *> &sites) {
	assert(sites.size() == LocalEne::fragmentlength);
	if (!defined(sites))
		return 0.0;

	std::string motif;
	std::string omigaseq;
	std::vector<double> tv;
//	std::cout <<sites[1]->resseq<<std::endl;
	for (int i = 1; i <= motiflength_; ++i) {
		const BackBoneSite &bs = *(sites[i]);
		if (bs.resname == "GLY")
			motif = motif + "G";
		else if (bs.resname == "PRO")
			motif = motif + "P";
		else
			motif = motif + "X";
		if (bs.omiga() < -90 || bs.omiga() > 90) {
			omigaseq += "c";
		} else {
			omigaseq += "t";
		}

		tv.push_back(bs.phi());
		tv.push_back(bs.psi());
	}
	std::string motifname = motif + omigaseq;
	return TorsionVectorScorer::getinstance().score(motifname, tv);
}
TetraSefEne::TetraSefEne() {
	std::string datapath = NSPdataio::datapath();
	RMinSEF::getinstance(datapath + "sstetra.dat", datapath + "coiltetra.dat");
	minsep_ = 7;
}
double TetraSefEne::energy(const std::vector<const BackBoneSite *> &frag1,
		const std::vector<const BackBoneSite *> &frag2, int sep) {
	if (sep < minsep_)
		return 0.0;
	int posim = LocalEne::fragmentlength / 2;
	const BackBoneSite &s1 = *frag1[posim];
	const BackBoneSite &s2 = *frag2[posim];
	return RMinSEF::getinstance().twobody(s1, s2);
}
double NSPproteinrep::calctotalenergy(MainChain &mc,
		const std::vector<SiteEne *> & terms,
		std::map<std::string, double> & energies, double egiveup) {
	for (auto iter = terms.begin(); iter != terms.end(); ++iter) {
		if (energies.find((*iter)->termname()) == energies.end())
			energies.insert(std::make_pair((*iter)->termname(), 0.0));
	}
	int posi = 0;
	for (auto &bs : mc) {
		if (mc.maskedormissing(posi++) != 0)
			continue;
		for (auto iter = terms.begin(); iter != terms.end(); ++iter) {
			energies[(*iter)->termname()] += (*iter)->energy(bs);
			if (energies[(*iter)->termname()] >= egiveup)
				return egiveup + 100.0;
		}
	}
	double res = 0.0;
	for (auto iter = terms.begin(); iter != terms.end(); ++iter) {
		energies[(*iter)->termname()] *= (*iter)->weight();
		res += energies[(*iter)->termname()];
	}
	return res;
}
double NSPproteinrep::calcloopenergy(MainChain &mc, int loopstart,
		int looplength, const std::vector<BackBoneSite> &newloop,
		const std::vector<SiteEne *> & terms,
		std::map<std::string, double> & energies, double egiveup) {
	for (auto iter = terms.begin(); iter != terms.end(); ++iter) {
		if (energies.find((*iter)->termname()) == energies.end())
			energies.insert(std::make_pair((*iter)->termname(), 0.0));
	}
	for (int i = 0; i < looplength; ++i) {
		const BackBoneSite *bsp;
		if (newloop.empty()) {
			bsp = &(mc.at(i + loopstart));
		} else {
			bsp = &(newloop[i]);
		}
//		if (mc.maskedormissing() != 0)
//			continue;
		for (auto iter = terms.begin(); iter != terms.end(); ++iter) {
			energies[(*iter)->termname()] += (*iter)->energy(*bsp);
			if (energies[(*iter)->termname()] >= egiveup)
				return egiveup + 100.0;
		}
	}
	double res = 0.0;
	for (auto iter = terms.begin(); iter != terms.end(); ++iter) {
		energies[(*iter)->termname()] *= (*iter)->weight();
		res += energies[(*iter)->termname()];
	}
	return res;
}
double NSPproteinrep::calctotalenergy(MainChain &mc,
		const std::vector<SitePairEne *> & terms,
		std::map<std::string, double> & energies, double egiveup) {
	for (auto iter = terms.begin(); iter != terms.end(); ++iter) {
		if (energies.find((*iter)->termname()) == energies.end())
			energies.insert(std::make_pair((*iter)->termname(), 0.0));
	}
	int posii = 0;
	for (auto iter1 = mc.begin(); iter1 != mc.end() - 1; ++iter1) {
		if (mc.maskedormissing(posii++) != 0)
			continue;
		auto &bs1 = *iter1;
		for (auto iter2 = iter1 + 1; iter2 != mc.end(); ++iter2) {
			auto &bs2 = *iter2;
			int sep = iter2 - iter1;
			if (mc.maskedormissing(posii + sep - 1) != 0)
				continue;
			for (auto iter = terms.begin(); iter != terms.end(); ++iter) {
//				double e=(*iter)->energy(bs1, bs2,sep);
//				if(e >10.0) {
//								std::cout << bs1.toString() << bs2.toString();
//				}
				energies[(*iter)->termname()] += (*iter)->energy(bs1, bs2, sep);
				if (energies[(*iter)->termname()] >= egiveup)
					return egiveup + 100.0;
			}
		}
	}
	double res = 0.0;
	for (auto iter = terms.begin(); iter != terms.end(); ++iter) {
		energies[(*iter)->termname()] *= (*iter)->weight();
		res += energies[(*iter)->termname()];
	}
	return res;
}
double NSPproteinrep::calcloopenergy(MainChain &mc, int loopstart,
		int looplength, const std::vector<BackBoneSite> &newloop,
		const std::vector<SitePairEne *> & terms,
		std::map<std::string, double> & energies, double egiveup) {
	for (auto iter = terms.begin(); iter != terms.end(); ++iter) {
		if (energies.find((*iter)->termname()) == energies.end())
			energies.insert(std::make_pair((*iter)->termname(), 0.0));
	}
	for (int i = 0; i < mc.size(); ++i) {
		bool usenew = (i >= loopstart && i < loopstart + looplength
				&& !newloop.empty());
		const BackBoneSite *bsp1;
		if (usenew)
			bsp1 = &newloop[i - loopstart];
		else {
			if (mc.maskedormissing(i) != 0)
				continue;
			bsp1 = &(mc[i]);
		}
		const BackBoneSite & bs1 = *bsp1;
		for (int j = 0; j < looplength; ++j) {
			if (i >= loopstart && i < loopstart + looplength) {
				if (i >= j + loopstart)
					continue;
			}
			const BackBoneSite *bsp;
			if (newloop.empty()) {
				bsp = &(mc.at(j + loopstart));
			} else {
				bsp = &(newloop[j]);
			}
			auto &bs2 = *bsp;
			int sep = j + loopstart - i;
			for (auto iter = terms.begin(); iter != terms.end(); ++iter) {
				if (sep > 0) {
					energies[(*iter)->termname()] += (*iter)->energy(bs1, bs2,
							sep);
				} else {
					energies[(*iter)->termname()] += (*iter)->energy(bs2, bs1,
							-sep);
				}
				if (energies[(*iter)->termname()] >= egiveup) {
//					std::cout <<"i-j,sep: " <<i <<"\t" <<j <<"\t"<<sep <<std::endl;
//					std::cout <<bs1.toString()<<bs2.toString();
					return egiveup + 100.0;
				}
			}
		}
	}
	double res = 0.0;
	for (auto iter = terms.begin(); iter != terms.end(); ++iter) {
		energies[(*iter)->termname()] *= (*iter)->weight();
		res += energies[(*iter)->termname()];
	}
	return res;
}
double NSPproteinrep::calctotalenergy(MainChain &mc,
		const std::vector<LocalEne *> & localterms,
		const std::vector<PackingEne *> & packingterms,
		std::map<std::string, double> & energies, double egiveup) {
	for (auto iter = localterms.begin(); iter != localterms.end(); ++iter) {
		if (energies.find((*iter)->termname()) == energies.end())
			energies.insert(std::make_pair((*iter)->termname(), 0.0));
	}
	for (auto iter = packingterms.begin(); iter != packingterms.end(); ++iter) {
		if (energies.find((*iter)->termname()) == energies.end())
			energies.insert(std::make_pair((*iter)->termname(), 0.0));
	}
	int posim = LocalEne::fragmentlength / 2;
	if (2 * posim == LocalEne::fragmentlength)
		posim -= 1;
	for (int i = 0; i < mc.size(); ++i) {
		if (mc.maskedormissing(i) != 0)
			continue;
		std::vector<const BackBoneSite *> fragment1;
		int posi0 = i - posim;
		for (int l = 0; l < LocalEne::fragmentlength; ++l) {
			int posi = posi0 + l;
			if (mc.missing(posi))
				fragment1.push_back(nullptr);
			else
				fragment1.push_back(&(mc.at(posi)));
		}
		if (!localterms.empty()) {
			for (auto iter = localterms.begin(); iter != localterms.end();
					++iter) {
				if (!(*iter)->defined(fragment1))
					continue;
				energies[(*iter)->termname()] += (*iter)->energy(fragment1)
						/ (*iter)->totalredundancy();
				if (energies[(*iter)->termname()] >= egiveup)
					return egiveup + 100.0;
			}
		}
		if (packingterms.empty())
			continue;
		int seqi = mc.at(i).resseq;
		for (int j = i + 1; j < mc.size(); ++j) {
			if (mc.maskedormissing(j) != 0)
				continue;
			std::vector<const BackBoneSite *> fragment2;
			int posij0 = j - posim;
			for (int l = 0; l < LocalEne::fragmentlength; ++l) {
				int posi = posij0 + l;
				if (mc.missing(posi))
					fragment2.push_back(nullptr);
				else
					fragment2.push_back(&(mc.at(posi)));
			}
			int sep = mc.at(j).resseq - seqi;
			for (auto iter = packingterms.begin(); iter != packingterms.end();
					++iter) {
				energies[(*iter)->termname()] += (*iter)->energy(fragment1,
						fragment2, sep);
				if (energies[(*iter)->termname()] >= egiveup)
					return egiveup + 100.0;
			}
		}
	}
	double res = 0.0;
	for (auto iter = localterms.begin(); iter != localterms.end(); ++iter) {
		energies[(*iter)->termname()] *= (*iter)->weight();
		res += energies[(*iter)->termname()];
	}
	for (auto iter = packingterms.begin(); iter != packingterms.end(); ++iter) {
		energies[(*iter)->termname()] *= (*iter)->weight();
		res += energies[(*iter)->termname()];
	}
	return res;
}
double NSPproteinrep::calcloopenergy(MainChain &mc, int loopstart,
		int looplength, const std::vector<BackBoneSite> &newloop,
		const std::vector<LocalEne *> & localterms,
		const std::vector<PackingEne *> & packingterms,
		std::map<std::string, double> & energies, double egiveup) {
	int posim = LocalEne::fragmentlength / 2;
	int coverstart = loopstart - posim;
	if (coverstart < 0)
		coverstart = 0;
	int coverend = loopstart + looplength + posim;
	if (coverend > mc.size())
		coverend = mc.size();
	int coverlength = coverend - coverstart;
	for (auto iter = localterms.begin(); iter != localterms.end(); ++iter) {
		if (energies.find((*iter)->termname()) == energies.end())
			energies.insert(std::make_pair((*iter)->termname(), 0.0));
	}
	for (auto iter = packingterms.begin(); iter != packingterms.end(); ++iter) {
		if (energies.find((*iter)->termname()) == energies.end())
			energies.insert(std::make_pair((*iter)->termname(), 0.0));
	}

	if (2 * posim == LocalEne::fragmentlength)
		posim -= 1;
//site_ene analysis
//	std::vector<double> locale_ene_site(coverend-coverstart,0.0);
//	std::vector<double> packing_ene_site(coverend-coverstart,0.0);

	for (int i = coverstart; i < coverend; ++i) {
		if (i < loopstart || i >= loopstart + looplength)
			if (mc.maskedormissing(i))
				continue;
		std::vector<const BackBoneSite *> fragment1;
		int posi0 = i - posim;
		for (int l = 0; l < LocalEne::fragmentlength; ++l) {
			int posi = posi0 + l;
			bool usenew = posi >= loopstart && posi < loopstart + looplength
					&& !newloop.empty();
			if (usenew) {
				fragment1.push_back(&newloop[posi - loopstart]);
			} else {
				if (mc.missing(posi))
					fragment1.push_back(nullptr);
				else
					fragment1.push_back(&(mc.at(posi)));
			}
		}
		if (!localterms.empty()) {
//			std::cout <<"loop "<< i <<std::endl;
			for (auto iter = localterms.begin(); iter != localterms.end();
					++iter) {
				if (!(*iter)->defined(fragment1))
					continue;
				if (i < loopstart - (*iter)->leftextension())
					continue;
				if (i >= loopstart + looplength + (*iter)->rightextension())
					continue;
				energies[(*iter)->termname()] += (*iter)->energy(fragment1)
						/ (*iter)->totalredundancy();
				if (energies[(*iter)->termname()] >= egiveup)
					return egiveup + 100.0;
			}
		}
		if (packingterms.empty())
			continue;
		if (i < loopstart || i >= loopstart + looplength)
			continue;
		int seqi = mc.at(i).resseq;
		for (int j = 0; j < mc.size(); ++j) {
			if (j >= loopstart && j < loopstart + looplength) {
				if (j > i)
					continue;
			} else if (mc.maskedormissing(j) != 0)
				continue;
			std::vector<const BackBoneSite *> fragment2;
			int posij0 = j - posim;
			for (int l = 0; l < LocalEne::fragmentlength; ++l) {
				int posi = posij0 + l;
				bool usenew = posi >= loopstart && posi < loopstart + looplength
						&& !newloop.empty();
				if (usenew) {
					fragment2.push_back(&newloop[posi - loopstart]);
				} else {
					if (mc.missing(posi))
						fragment2.push_back(nullptr);
					else
						fragment2.push_back(&(mc.at(posi)));
				}
			}
			int sep = mc.at(j).resseq - seqi;
			for (auto iter = packingterms.begin(); iter != packingterms.end();
					++iter) {
				if (sep > 0) {
					energies[(*iter)->termname()] += (*iter)->energy(fragment1,
							fragment2, sep);
				} else {
					energies[(*iter)->termname()] += (*iter)->energy(fragment2,
							fragment1, -sep);
				}
				if (energies[(*iter)->termname()] >= egiveup)
					return egiveup + 100.0;
			}
		}
	}
	double res = 0.0;
	for (auto iter = localterms.begin(); iter != localterms.end(); ++iter) {
		energies[(*iter)->termname()] *= (*iter)->weight();
		res += energies[(*iter)->termname()];
	}
	for (auto iter = packingterms.begin(); iter != packingterms.end(); ++iter) {
		energies[(*iter)->termname()] *= (*iter)->weight();
		res += energies[(*iter)->termname()];
	}
	return res;
}
double NSPproteinrep::calctotalenergy(MainChain &mc, EnergyTerms &terms,
		std::map<std::string, double> & energies, double egiveup) {
	double ene = 0.0;
	if (!terms.siteene_terms.empty()) {
		ene += calctotalenergy(mc, terms.siteene_terms, energies, egiveup);
		if (ene >= egiveup)
			return ene;
	}
	if (!terms.sitepairene_terms.empty()) {
		ene += calctotalenergy(mc, terms.sitepairene_terms, energies, egiveup);
		if (ene >= egiveup)
			return ene;
	}
	if (!terms.localene_terms.empty() || !terms.packingene_terms.empty()) {
		ene += calctotalenergy(mc, terms.localene_terms, terms.packingene_terms,
				energies, egiveup);
	}
	return ene;
}
double NSPproteinrep::calcloopenergy(MainChain &mc, int loopstart,
		int looplength, const std::vector<BackBoneSite> &newloop,
		EnergyTerms &terms, std::map<std::string, double> & energies,
		double egiveup) {
	assert(loopstart + looplength <= mc.size());
	double ene = 0.0;
	if (!terms.siteene_terms.empty()) {
		ene += calcloopenergy(mc, loopstart, looplength, newloop,
				terms.siteene_terms, energies, egiveup);
		if (ene >= egiveup)
			return ene;
	}
	if (!terms.sitepairene_terms.empty()) {
		ene += calcloopenergy(mc, loopstart, looplength, newloop,
				terms.sitepairene_terms, energies, egiveup);
		if (ene >= egiveup)
			return ene;
	}
	if (!terms.localene_terms.empty() || !terms.packingene_terms.empty()) {
		ene += calcloopenergy(mc, loopstart, looplength, newloop,
				terms.localene_terms, terms.packingene_terms, energies,
				egiveup);
	}
	return ene;
}
double NSPproteinrep::looploopenergy(const std::vector<BackBoneSite> &loop1,
		const std::vector<BackBoneSite> &loop2, SitePairEne * eterm,
		double egiveup) {
	double ene = 0.0;
	for (auto iter1 = loop1.begin(); iter1 != loop1.end(); ++iter1) {
		auto &bs1 = *iter1;
		for (auto iter2 = loop2.begin(); iter2 != loop2.end(); ++iter2) {
			auto &bs2 = *iter2;
			ene += eterm->energy(bs1, bs2, 1000);
			if (ene >= egiveup)
				return egiveup + 100.0;
		}
	}
	return ene;
}
void NSPproteinrep::decomposeenergy(MainChain &mc,const EnergyTerms & eterms,std::ostream &os) {
	const std::vector<LocalEne *> & localterms=eterms.localene_terms;
	const std::vector<PackingEne *> &packingterms=eterms.packingene_terms;
	int posim = LocalEne::fragmentlength / 2;
	if (2 * posim == LocalEne::fragmentlength)
		posim -= 1;
	for (int i = 0; i < mc.size(); ++i) {
		if (mc.maskedormissing(i) != 0)
			continue;
		std::vector<const BackBoneSite *> fragment1;
		int posi0 = i - posim;
		for (int l = 0; l < LocalEne::fragmentlength; ++l) {
			int posi = posi0 + l;
			if (mc.missing(posi))
				fragment1.push_back(nullptr);
			else
				fragment1.push_back(&(mc.at(posi)));
		}
		if (!localterms.empty()) {
			for (auto iter = localterms.begin(); iter != localterms.end();
					++iter) {
				if (!(*iter)->defined(fragment1))
					continue;
				double etmp=(*iter)->energy(fragment1)
								/ (*iter)->totalredundancy();
				os << i+1 <<"\t"<<(*iter)->termname()<<" " << etmp <<std::endl;
			}
		}
		if (packingterms.empty())
			continue;
		int seqi = mc.at(i).resseq;
		for (int j = 0 ; j < mc.size(); ++j) {
			if (mc.maskedormissing(j) != 0)
				continue;
			if (j==i) continue;
			std::vector<const BackBoneSite *> fragment2;
			int posij0 = j - posim;
			for (int l = 0; l < LocalEne::fragmentlength; ++l) {
				int posi = posij0 + l;
				if (mc.missing(posi))
					fragment2.push_back(nullptr);
				else
					fragment2.push_back(&(mc.at(posi)));
			}
			int sep = mc.at(j).resseq - seqi;
			for (auto iter = packingterms.begin(); iter != packingterms.end();
					++iter) {
				double etmp;
				if(sep > 0)
					etmp=(*iter)->energy(fragment1,
						fragment2, sep);
				else
					etmp=(*iter)->energy(fragment2,
											fragment1, -sep);
				os <<i+1 <<"-"<<j+1 <<"\t" <<(*iter)->termname()<<" "<<etmp<<std::endl;
			}
		}
	}
}
