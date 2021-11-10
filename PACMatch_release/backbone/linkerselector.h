/*
 * linkerselector.h
 *
 *  Created on: 2016年12月3日
 *      Author: hyliu
 */

#ifndef BACKBONE_LINKERSELECTOR_H_
#define BACKBONE_LINKERSELECTOR_H_
#include "backbone/backbonesite.h"
#include <vector>
#include <map>
namespace NSPproteinrep {
template<typename INTERACTION>
class LinkerSelector {
public:
	typedef std::vector<BackBoneSite> Segment;
	struct ScoredLinker {
		ScoredLinker(Segment * sg, double sc) :
				linkersegment(sg), linkerscore(sc) {
			;
		}
		Segment *linkersegment;
		double linkerscore;
	};

	LinkerSelector(std::vector<Segment> *sssegments) :
			SSsegments_(sssegments) {
		;
	}

	std::vector<unsigned int> & segmentseq() {
		return segmentseq_;
	}
	const std::vector<unsigned int> & segmentseq() const { return segmentseq_;}
	Segment & getSSsegment(unsigned int i) const {return SSsegments_->at(segmentseq_[i]);}
	void addLinker(unsigned int seg1, unsigned int seg2,
			Segment *linkersegment) {
		std::pair<unsigned int, unsigned int> key(seg1, seg2);
		if (linkers_.find(key) == linkers_.end())
			linkers_.insert(std::make_pair(key, std::vector<ScoredLinker>()));
		double score = INTERACTION::interaction(*linkersegment);
		for (unsigned int s = 0; s < SSsegments_->size(); ++s) {
			if (s == seg1)
				score += INTERACTION::interaction(SSsegments_->at(s), *linkersegment, true);
			else if (s == seg2)
				score += INTERACTION::interaction(*linkersegment, SSsegments_->at(s), true);
			else
				score += INTERACTION::interaction(SSsegments_->at(s), *linkersegment);
		}
		if (score < INTERACTION::SCORECUT)
			linkers_.at(key).push_back(ScoredLinker(linkersegment, score));
	}
	double operator()(const std::vector<Segment *> & prevlinkers,
			std::vector<std::pair<Segment *, double>> *nextlinkers) {
		unsigned int nstep = prevlinkers.size();
		if (segmentseq_.size() <= nstep + 1)
			return 0.0;
		unsigned int seg1 = segmentseq_[nstep];
		unsigned int seg2 = segmentseq_[nstep + 1];
		std::pair<unsigned int, unsigned int> key(seg1, seg2);
		auto it = linkers_.find(key);
		if (it == linkers_.end())
			return 0.0;
		double scoretot = 0.0;
		for (auto & sl : it->second) {
			double score = sl.linkerscore;
			for (auto & pl : prevlinkers) {
				score += INTERACTION::interaction(*pl, *(sl.linkersegment));
			}
			if (score < INTERACTION::SCORECUT) {
				score = exp(-score);
				nextlinkers->push_back(std::make_pair(sl.linkersegment, score));
				scoretot += score;
			}
		}
		return scoretot;
	}
private:
	std::vector<Segment> *SSsegments_;
	std::vector<unsigned int> segmentseq_;
	std::map<std::pair<unsigned int, unsigned int>, std::vector<ScoredLinker>> linkers_;
};

}

#endif /* BACKBONE_LINKERSELECTOR_H_ */
