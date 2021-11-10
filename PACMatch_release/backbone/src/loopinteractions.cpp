/*
 * linkerinteractions.cpp
 *
 *  Created on: 2016年12月5日
 *      Author: hyliu
 */

#include "backbone/loopinteractions.h"
#include "backbone/rminsef.h"
#include "backbone/torsionvectorscorer.h"
using namespace NSPproteinrep;

double LoopInteractions::SCORECUT { 1.e10 };
LoopInteractions::Switches LoopInteractions::switches { true, true, true,true};
LoopInteractions::Switches LoopInteractions::switchstate_steric_torsion { true, false,false,true};
LoopInteractions::Switches LoopInteractions::switchstate_torsion_motif_tetra {false,true,true,true};
typedef std::vector<BackBoneSite> Segment;
double LoopInteractions::steric_interaction(const Segment & seg) {
	double score = 0.0;
	for (auto & s : seg) {
		score += SiteInteractions::steric_interaction(s);
	}
	for (unsigned int i = 0; i < seg.size() - 1; ++i) {
		for (unsigned int j = i + 1; j < seg.size(); ++j) {
			score += SiteInteractions::steric_interaction(seg[i], seg[j],
					j - i);
			if (score > SCORECUT)
				return score;
		}
	}
	return score;
}
double LoopInteractions::steric_interaction(const Segment & seg1,
		const Segment & seg2, bool successive) {
	double score = 0.0;
	unsigned int sep1 = seg1.size();
	unsigned int posi = 0;
	for (auto & s1 : seg1) {
		unsigned int sep2 = 0;
		for (auto & s2 : seg2) {
			unsigned int sep = sep1 + sep2;
			if (successive)
				score += SiteInteractions::steric_interaction(s1, s2, sep);
			else
				score += SiteInteractions::steric_interaction(s1, s2);
			if (score > LoopInteractions::SCORECUT)
				return score;
			++sep2;
		}
		--sep1;
	}
	return score;
}

double LoopInteractions::sef_interaction(const Segment & seg1,
		const Segment & seg2, bool successive) {
	double score = 0.0;
	unsigned int sep1 = seg1.size();
	unsigned int posi = 0;
	RMinSEF &rminsef = RMinSEF::getinstance();
	for (auto & s1 : seg1) {
		unsigned int sep2 = 0;
		for (auto & s2 : seg2) {
			unsigned int sep = sep1 + sep2;
			if (sep > 6 || !successive)
				score += rminsef.twobody(s1, s2);
			++sep2;
		}
		--sep1;
	}
	return score;
}
double SiteInteractions::steric_interaction(const BackBoneSite & s) {
	return 0.0;
}

double SiteInteractions::steric_interaction(const BackBoneSite &s1,
		const BackBoneSite &s2) {
	if (atomsclashed(s1, s2)) {
//		std::cout <<"Clashed: " <<s1.resid<<s1.resname <<" " <<s2.resid<<s2.resname<<std::endl;
		return 1.e20;
	}
	return 0.0;
}

double SiteInteractions::steric_interaction(const BackBoneSite &s1,
		const BackBoneSite & s2, unsigned int sep) {
	if (sep != 1)
		return steric_interaction(s1, s2);
	return 0.0;
}
double NSPproteinrep::loop_context_steric_interaction(
		std::vector<std::vector<BackBoneSite>> *context,
		unsigned int headsegment, unsigned int tailsegment,
		const std::vector<BackBoneSite> & loop) {
	double score = LoopInteractions::steric_interaction(loop);
	for (unsigned int i = 0; i < context->size(); ++i) {
		if (i == headsegment) {
			std::vector<BackBoneSite> seg1;
			unsigned int size = (*context)[i].size() - 1;
			seg1.resize(size, BackBoneSite());
			std::copy((*context)[i].begin(), (*context)[i].begin() + size,
					seg1.begin());
			score += LoopInteractions::steric_interaction(seg1, loop, true);
			if (score >= LoopInteractions::SCORECUT)
				return score;
		} else if (i == tailsegment) {
			std::vector<BackBoneSite> seg2;
			unsigned int size = (*context)[i].size() - 1;
			seg2.resize(size, BackBoneSite());
			std::copy((*context)[i].begin() + 1, (*context)[i].end(),
					seg2.begin());
			score += LoopInteractions::steric_interaction(loop, seg2, true);
			if (score >= LoopInteractions::SCORECUT)
				return score;
		} else {
			score += LoopInteractions::steric_interaction((*context)[i], loop);
			if (score > LoopInteractions::SCORECUT)
				return score;
		}
	}
	return score;
}

double NSPproteinrep::loop_context_interaction(
		std::vector<std::vector<BackBoneSite>> *context,
		unsigned int headsegment, unsigned int tailsegment,
		const std::vector<BackBoneSite> & loop) {
	static bool initialized { false };
	if (!initialized) {
		if (LoopInteractions::switches.tetrasef_on)
			RMinSEF &rminsef = RMinSEF::getinstance("sstetra.dat",
					"coiltetra.dat");
		if (LoopInteractions::switches.torsionmotif_on) {
			std::vector<BackBoneSite> tmpsites;
			readbackbonesites("tmplatesites.dat", tmpsites);
			TorsionVectorScorer &tvscorer = TorsionVectorScorer::getinstance(
					&tmpsites);
		}
		initialized = true;
	}
	double score = 0.0;
	if (LoopInteractions::switches.steric_on)
		score += loop_context_steric_interaction(context, headsegment,
				tailsegment, loop);
	if (score > LoopInteractions::SCORECUT)
		return score;
	if (LoopInteractions::switches.torsionmotif_on) {
		std::vector<BackBoneSite> chain;
		int chainlength = (*context)[headsegment].size()
				+ (*context)[tailsegment].size() + loop.size() - 2;
		chain.resize(chainlength, BackBoneSite());
		int l1 = (*context)[headsegment].size() - 1;
		int l2 = l1 + loop.size();
		std::copy((*context)[headsegment].begin(),
				(*context)[headsegment].begin() + l1, chain.begin());
		std::copy(loop.begin(), loop.end(), chain.begin() + l1);
		std::copy((*context)[tailsegment].begin() + 1,
				(*context)[tailsegment].end(), chain.begin() + l2);
		TorsionVectorScorer &motifscorer = TorsionVectorScorer::getinstance();
		int scorebegin = l1 - motifscorer.length() + 1;
		if (scorebegin < 0)
			scorebegin = 0;
		int scoreend = l2 + motifscorer.length() - 1;
		if (scoreend > chain.size())
			scoreend = chain.size();
		score += motifscorer.scorerange(chain.begin() + scorebegin,
				chain.begin() + scoreend);
	}
	if (LoopInteractions::switches.tetrasef_on) {
		for (unsigned int i = 0; i < context->size(); ++i) {
			if (i == headsegment) {
				std::vector<BackBoneSite> seg1;
				unsigned int size = (*context)[i].size() - 1;
				seg1.resize(size, BackBoneSite());
				std::copy((*context)[i].begin(), (*context)[i].begin() + size,
						seg1.begin());
				score += LoopInteractions::sef_interaction(seg1, loop, true);
			} else if (i == tailsegment) {
				std::vector<BackBoneSite> seg2;
				unsigned int size = (*context)[i].size() - 1;
				seg2.resize(size, BackBoneSite());
				std::copy((*context)[i].begin() + 1, (*context)[i].end(),
						seg2.begin());
				score += LoopInteractions::sef_interaction(loop, seg2, true);
				if (score >= LoopInteractions::SCORECUT)
					return score;
			} else {
				score += LoopInteractions::sef_interaction((*context)[i], loop);
			}
		}
	}
	return score;

}


double NSPproteinrep::loop_loop_sef(std::vector<BackBoneSite>::const_iterator & iter_a,
		int posi_a, int length_a,
		std::vector<BackBoneSite>::const_iterator & iter_b, int posi_b,int length_b){
	double score = 0.0;
	RMinSEF &rminsef = RMinSEF::getinstance();
	for (int i=0; i<length_a;++i) {
		auto bsa=iter_a+i;
		for(int j=0; j<length_b;++j) {
			auto bsb=iter_b+j;
			int sep=(posi_b+j)-(posi_a+i);
			if(0<sep<=6) {
				score +=rminsef.twobody(*bsa, *bsb);
			} else if (-6<=sep<0) {
				score +=rminsef.twobody(*bsb,*bsa);
			}
		}
	}
	return score;
}
