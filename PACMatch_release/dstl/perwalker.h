/*
 * perwalker.h
 *
 *  Created on: 2016年12月2日
 *      Author: hyliu
 */

#ifndef DSTL_PERWALKER_H_
#define DSTL_PERWALKER_H_
#include "dstl/randomengine.h"
#include <vector>
namespace NSPdstl{
struct PerParameters{
	PerParameters(unsigned int l):length(l),branch(1.0),enrich(1.0),attrition(1.0){;}
	std::vector<std::pair<double,double>>  logweightlimits;
	double branch;
	double enrich;
	double attrition;
	unsigned int length;
};

/*
 * prunning-enrichment walker
 */
template <typename BEAD>
class PerWalker{
public:
	PerWalker(const PerParameters & p): par_(p) {;}
	PerWalker(unsigned int l):par_(l){;}
	template <typename NEXTBEADSELECTOR, typename ATLEAF>
	void permcstep(NEXTBEADSELECTOR & nextselector, ATLEAF & atleaf) {
		unsigned int nstep=beads_.size();
		if(nstep == par_.length) {
			atleaf(*this);
			return;
		}
		std::vector<std::pair<BEAD,double>> nextcandidates;
		double ztot=nextselector(beads_,&nextcandidates);
		if(nextcandidates.size() == 0) {
			return;
		}
		double newlw=0.0;
		double enrich=1.0;
		if(nstep >0) {
			newlw=logweights_.back()+logweightshifts_.back()+ztot;
			enrich=par_.branch;
			if(newlw >par_.logweightlimits[nstep].second) {
				enrich *=par_.enrich;
			} else if(newlw <par_.logweightlimits[nstep].first){
				enrich *=par_.attrition;
			}
		}
		RandomEngine<> & re= RandomEngine<>::getinstance();
		re.setrealrng(0.0,1.0);
		while (re.randomreal() < enrich ) {
			double random=re.randomreal()*ztot;
			double zpart=0.0;
			for (auto & b:nextcandidates) {
				zpart += b.second;
				if(zpart >= random ) {
					moveforward(b,newlw,enrich);
					permcstep(nextselector,atleaf);
					rollback();
					break;
				}
			}
			enrich -=1.0;
		}
		return;
	}

	template <typename NEXTBEADSELECTOR, typename ATLEAF>
	void enumstep(NEXTBEADSELECTOR & nextselector, ATLEAF & atleaf) {
		if(beads_.size() == par_.length) {
			atleaf(*this);
			return;
		}
		std::vector<std::pair<BEAD,double>> nextcandidates;
		nextselector(beads_,&nextcandidates);
		if(nextcandidates.size() == 0) {
			return;
		}
		for (auto & b:nextcandidates) {
				if(b.second <=0.0) continue;
				moveforward(b,0.0,1.0);
				enumstep(nextselector,atleaf);
				rollback();
		}
		return;
	}
	void moveforward(std::pair<BEAD,double> & b, double newlw,double enrich) {
		beads_.push_back(b.first);
		enrich=-log(enrich);
		if(!logweights_.empty()) {
			newlw += logweights_.back();
			enrich += logweightshifts_.back();
		}
		logweights_.push_back(newlw);
		logweightshifts_.push_back(enrich);
		bweights_.push_back(b.second);
	}
	void rollback() {
		beads_.pop_back();
		logweights_.pop_back();
		logweightshifts_.pop_back();
		bweights_.pop_back();
	}
	std::vector<BEAD> & beads() {return beads_;}
	const std::vector<BEAD> & beads() const {return beads_;}
private:
	std::vector<BEAD> beads_;
	std::vector<double> logweights_;
	std::vector<double> logweightshifts_;
	std::vector<double> bweights_;
	PerParameters par_;
};

}



#endif /* DSTL_PERWALKER_H_ */
