/*
 * testtvscore.cpp
 *
 *  Created on: 2016年12月28日
 *      Author: hyliu
 */

#include "backbone/torsionvectorscorer.h"
#include "dstl/randomengine.h"
#include "pdbstatistics/phipsidistr.h"
using namespace NSPproteinrep;
std::string randomglypro() {
	NSPdstl::RandomEngine<> &rneg = NSPdstl::RandomEngine<>::getinstance();
	const double PGLY { 0.024 };
	const double PPRO { 0.049 };
	double p = rneg.realrng(0.0, 1.0)();
	if (p < PGLY) {
		return "GLY";
	} else if (p < PGLY + PPRO)
		return "PRO";
	return "RES";
}
void randomtv(int length, std::string *motifname, std::vector<double> *tv) {
	std::vector<BackBoneSite> sites(length + 1, BackBoneSite());
	for (auto &bs : sites) {
		bs.resname = randomglypro();
		bs.data_[BackBoneSite::OMIGA] = 180.0;
	}
	NSPdstl::RandomEngine<> &rneg = NSPdstl::RandomEngine<>::getinstance();
	for (int i = 0; i < length; ++i) {
		if (sites[i + 1].resname == "PRO") {
			double u = rneg.realrng(0.0, 1.0)();
			if (u < 0.1)
				sites[i].data_[BackBoneSite::OMIGA] = 0.0;
		}
	}
	bool iscis = false;
	for (int i = 0; i < length; ++i) {
		const NSPpdbstatistics::PhiPsiDistr *distr;
		distr = &(NSPpdbstatistics::PhiPsiDistr::coildistr());
		if (sites[i].resname == "GLY") {
			distr = &(NSPpdbstatistics::PhiPsiDistr::glydistr());
		} else if (sites[i].resname == "PRO") {
			if (iscis)
				distr = &(NSPpdbstatistics::PhiPsiDistr::cisprodistr());
			else
				distr = &(NSPpdbstatistics::PhiPsiDistr::transprodistr());
		} else {
			if (sites[i + 1].resname == "PRO")
				distr = &(NSPpdbstatistics::PhiPsiDistr::preprodistr());
		}
		distr->randomphipsi(rneg.realrng(),
				&(sites[i].data_[BackBoneSite::PHI]),
				&(sites[i].data_[BackBoneSite::PSI]));
		if (sites[i].omiga() < 90.0)
			iscis = true;
		else
			iscis = false;
	}
	*motifname = getmotifname(sites.begin(), length);
	tv->clear();
	for (int i = 0; i < length; ++i) {
		tv->push_back(sites[i].phi());
		tv->push_back(sites[i].psi());
	}
}
int main(int argc, char **argv) {
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
	TorsionVectorScorer tvscorer;
	int length = 4;
	tvscorer.init(&sites, length);

	std::vector<BackBoneSite> testsites;
	readbackbonesites(std::string(argv[2]), testsites);
	std::cout << "Finish reading testsites" << std::endl;
	double sumnt = 0.0;
	double sumnt2 = 0.0;
	double countnt = 0.0;
	int i = 0;
	while (countnt < 1000.0) {
		i +=10;
		std::vector<BackBoneSite>::iterator bs = testsites.begin() + i;
		if (!fragstartsite(bs, testsites.end(), length))
			continue;
		bool containcoil = false;
		for (int ii = 0; ii < length; ++ii) {
			auto it = bs + ii;
			if (it->sscodechar() == 'C') {
				containcoil = true;
				break;
			}
		}
		if (!containcoil)
			continue;
		std::string mtype=NSPproteinrep::getmotifname(bs,length);
		if (mtype.find('P') != std::string::npos //) continue;
		||		mtype.find('G') != std::string::npos) continue;
		double score = tvscorer.score(bs);
		sumnt += score;
		countnt += 1.0;
		sumnt2 += score * score;
		std::cout << "native " << i << "\t" <<mtype<<"\t" << score << std::endl;
	}
	double sumrnd = 0.0;
	double sumrnd2 = 0.0;
	double countrnd = 0.0;
	while (countrnd<1000.0) {
		std::string motifname;
		std::vector<double> tv;
		randomtv(length, &motifname, &tv);
		if(motifname.find('P') != std::string::npos || //) continue;
				motifname.find('G') != std::string::npos) continue;
		double score = tvscorer.score(motifname, tv);
		std::cout << "random " << i << "\t" <<motifname<<"\t"<< score << std::endl;
		sumrnd += score;
		sumrnd2 += score * score;
		countrnd += 1.0;
	}

	sumrnd /= countrnd;
	std::cout << "averagelogscorernd: " << sumrnd << "\t"
			<< sqrt(sumrnd2 / countrnd - (sumrnd * sumrnd)) << std::endl;
	sumnt /= countnt;
	std::cout << "averagelogscorent: " << sumnt<< "\t"
			<< sqrt(sumnt2 / countnt - (sumnt * sumnt)) << std::endl;
}
