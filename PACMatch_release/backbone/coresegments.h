/*
 * coresegments.h
 *
 *  Created on: 2017年5月11日
 *      Author: hyliu
 */

#ifndef BACKBONE_CORESEGMENTS_H_
#define BACKBONE_CORESEGMENTS_H_
#include "backbone/backbonesite.h"
#include "backbone/mainchain.h"
#include <iostream>
#include <set>
#include <map>
#include <memory>
namespace NSPproteinrep{
typedef std::vector<BackBoneSite> Segment;
typedef std::shared_ptr<Segment> SegmentPtr;
class CoreSegments: public std::vector<Segment>{
public:
	struct LoopSolution {
		SegmentPtr loop;
		int nextension{0};
		int cextension{0};
		double energy{0};
	};
	struct LinkedChain {
		std::vector<int> segments;
		std::vector<int> loopidx_;
		double energy{0.0};
	};
	typedef std::pair<int,int> GapIndex;
	typedef std::vector<LoopSolution> LoopSolutions;
	bool read(std::istream &is);
	int nflexible(int s) const {return n_flexible_[s];}
	int cflexible(int s) const {return c_flexible_[s];}
	int sizeofsegment(int s) const {return this->at(s).size();}
	int buildclosedloops(int s1,int s2,int minlength=1,int maxlength=10);
	int buildlinkedchains(const std::vector<int> & order);
	void extendchain(LinkedChain lc);
	int linkedchainlength(const LinkedChain & lc) const;
	int linkedchainlooplength(const LinkedChain &lc) const;
	double linkedchainaverageene(const LinkedChain &lc) const;
	std::string linkedchainname(const LinkedChain &lc) const;
	void buildlinkedmainchain(const LinkedChain &lc, MainChain *mc) const;
	std::vector<LinkedChain> & linkedchains() {return linkedchains_;}
	const std::vector<LinkedChain> & linkedchains() const {return linkedchains_;}
	const LoopSolutions * loopsbetween(int s1,int s2) const {
		GapIndex idx=std::make_pair(s1,s2);
		auto it=closedloops_.find(idx);
		if( it != closedloops_.end()) return &(it->second);
		else return nullptr;
	}
	void sortlinkedchains();
private:
	std::vector<int> n_flexible_;
	std::vector<int> c_flexible_;
	std::map<GapIndex,LoopSolutions> closedloops_;
	void buildmainchain(const std::vector<int> & order,
				const std::vector<int> & gaps, MainChain *mc) const;
	std::vector<LinkedChain> linkedchains_;
	std::set<std::string> builtchains_;
};
}


#endif /* BACKBONE_CORESEGMENTS_H_ */
