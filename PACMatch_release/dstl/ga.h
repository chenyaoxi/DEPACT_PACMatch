/*
 * ga.h
 *
 *  Created on: 2017年11月13日
 *      Author: hyliu
 */

#ifndef DSTL_GA_H_
#define DSTL_GA_H_
#include "dstl/sortindex.h"
#include "dstl/randomengine.h"
#include <vector>
#include <memory>
#include <cassert>
namespace NSPdstl {
template<typename SPECIES>
class GAPopulation {
public:
	typedef std::shared_ptr<SPECIES> IVPtr;  //point to an individual
	GAPopulation(int size, int subsize) :
			populationsize_(size), subpopsize_(subsize) {
		rng_ = std::shared_ptr<RandomEngine<>>(new RandomEngine<>);
		rng_->init();
	}
	void addindividual(IVPtr iv, double fitness) {
		individuals_.push_back(iv);
		fitness_.push_back(fitness);
	}
	const std::vector<IVPtr> & individuals() const {
		return individuals_;
	}
	const std::vector<double> & fitness() const {
		return fitness_;
	}
	void selection(int newsize) {
		if (newsize >= fitness_.size())
			return;
		sort12(fitness_, individuals_);
		fitness_.resize(newsize);
		individuals_.resize(newsize);
	}
	int ncenter() const {return centers_.size();}
	void selection_noise(const std::vector<int> & sizes,
			std::vector<double> &noisesds) {
		int newsize = 0;
		for (auto s : sizes)
			newsize += s;
		if (newsize >= fitness_.size())
			return;
		std::vector<int> indices;
		for (int i = 0; i < fitness_.size(); ++i)
			indices.push_back(i);
		int start = 0;
		for (int m = 0; m < sizes.size(); ++m) {
			std::vector<double> noisefit;
			std::vector<int> sortindex;
			for (int i = start; i < fitness_.size(); ++i) {
				noisefit.push_back(fitness_[indices[i]]);
				sortindex.push_back(indices[i]);
			}
			for (auto &f : noisefit)
				f += rng_->randomnormal(0.0, noisesds[m]);
			sort12(noisefit, sortindex);
			for (int i = start; i < fitness_.size(); ++i) {
				indices[i] = sortindex[i - start];
			}
			start += sizes[m];
		}
		std::vector<double> fitness_old = fitness_;
		std::vector<IVPtr> ivold = individuals_;
		fitness_.resize(newsize);
		individuals_.resize(newsize);
		for (int i = 0; i < newsize; ++i) {
			fitness_[i] = fitness_old[indices[i]];
			individuals_[i] = ivold[indices[i]];
		}
	}
	void selection(const std::vector<std::pair<int, double>> &sizenshifts) {
		sort12(fitness_, individuals_);
		double emin = fitness_[0];
		std::vector<int> ntotal;
		double elow = emin;
		for (int i = 0; i < sizenshifts.size(); ++i) {
			ntotal.push_back(0);
			int & nt = ntotal.back();
			double ehigh = emin + sizenshifts[i].second;
			for (auto f : fitness_)
				if (f > elow && f <= ehigh)
					++nt;
			elow = ehigh;
		}
		int start = 1;
		std::vector<bool> keep(fitness_.size(), true);
		for (int i = 0; i < sizenshifts.size(); ++i) {
			if (ntotal[i] <= sizenshifts[i].first) {
				start += ntotal[i];
				continue;
			}
			int nremove = ntotal[i] - sizenshifts[i].first;
			int nremoved = 0;
			double premove = nremove / ntotal[i];
			while (nremoved < nremove) {
				int iremove = start + rng_->intrng(0, ntotal[i] - 1)();
				if (keep[iremove]) {
					keep[iremove] = false;
					++nremoved;
				}
			}
			start += ntotal[i];
		}
		std::vector<double> fitness_old = fitness_;
		std::vector<IVPtr> ivold = individuals_;
		fitness_.clear();
		individuals_.clear();
		for (int i = 0; i < start; ++i) {
			if (keep[i]) {
				fitness_.push_back(fitness_old[i]);
				individuals_.push_back(ivold[i]);
			}
		}
	}
	template<typename DISTFUNC>
	int selection_distance(const DISTFUNC &distscore) {
//		if (fitness_.size() <= populationsize_)
//			return 0;
		sort12(fitness_, individuals_);
		double emin = fitness_[0];
		std::vector<double> sortingscore(fitness_.size(), 0.0);
		centers_.clear();
		std::vector<int> neighbors;
		centers_.push_back(std::make_pair(individuals_[0], fitness_[0]));
		neighbors.assign(fitness_.size(), 0);
		double scoremin = distscore(individuals_[0], individuals_[0]);
		sortingscore[0] = fitness_[0] + scoremin-10000;
		for (int i = 1; i < fitness_.size(); ++i) {
			double smax = -1.e10;
			bool asneighbor = false;
			for (int ic = 0; ic < centers_.size(); ++ic) {
				double sic = distscore(individuals_[i], centers_[ic].first);
				if (sic > scoremin) {
					++neighbors[ic];
					if (neighbors[ic] < subpopsize_)
						asneighbor = true;
				}
				if (sic > smax)
					smax = sic;
			}

			if (smax == scoremin) {
				centers_.push_back(
						std::make_pair(individuals_[i], fitness_[i]));
				smax -=10000; //centers are considered to keep first
			}
			if (asneighbor)
				smax = scoremin;
			sortingscore[i] = smax + fitness_[i];
		}
		std::vector<int> resortedindex = sortindex(sortingscore);
		std::vector<double> fitness_old = fitness_;
		std::vector<IVPtr> ivold = individuals_;
		if(fitness_.size()>populationsize_){
			fitness_.resize(populationsize_);
			individuals_.resize(populationsize_);
		}
		for (int i = 0; i < fitness_.size(); ++i) {
			fitness_[i] = fitness_old[resortedindex[i]];
			individuals_[i] = ivold[resortedindex[i]];
		}
		return centers_.size();
	}
	template<typename DISTFUNC>
	void StoreNShrink(const DISTFUNC &distscore) {
		assert(!centers_.empty());
		double scoremin = distscore(centers_[0].first, centers_[0].first);
		std::vector<bool> keep(centers_.size(), true);
		std::vector<bool> store(centers_.size(),false);
		for (int i = 0; i < centers_.size(); ++i) {
			for (int j = 0; j < visited_.size(); ++j) {
				if (centers_[i].first == visited_[j].first) {
					keep[i] = false;
					break;
				}
			}
		}
		std::vector<std::pair<int,int>> replace;
		for(int i=0;i<centers_.size();++i){
			if (!keep[i])
				continue;
			store[i]=true;
			for (int j = 0; j < visited_.size(); ++j) {
				double sic = distscore(centers_[i].first, visited_[j].first);
				if (sic > scoremin) {
					store[i]=false;
					if (centers_[i].second < visited_[j].second) {
						replace.push_back(std::make_pair(i,j));
						break;
					} else {
						keep[i]=false;
					}
				}
			}
		}
		for(auto &r:replace){
			visited_[r.second]=centers_[r.first];
		}
		fitness_.clear();
		individuals_.clear();
		int maxkeep=(int)(populationsize_/subpopsize_)*0.618;
		if(maxkeep<1) maxkeep=1;
		int nkeep=0;
		for(int i=0;i<centers_.size();++i){
			if(keep[i]){
				++nkeep;
				if(nkeep<=maxkeep){
					individuals_.push_back(centers_[i].first);
					fitness_.push_back(centers_[i].second);
				} else {
					break;
				}
			}
			if(store[i]){
				visited_.push_back(centers_[i]);
			}
		}
		assert(nkeep>0);
		auto f=[](const StoredIV &x,const StoredIV & y)->bool{return x.second < y.second;};
		std::sort(visited_.begin(),visited_.end(),f);
	}
	const std::vector<std::pair<IVPtr,double>> &visited()const {return visited_;}
private:
	typedef std::pair<IVPtr, double> StoredIV;
	std::vector<double> fitness_;
	std::vector<IVPtr> individuals_;
	int populationsize_;
	int subpopsize_;
	std::vector<StoredIV> centers_;
	std::vector<StoredIV> visited_;
	std::shared_ptr<NSPdstl::RandomEngine<>> rng_ { nullptr };
};

}

#endif /* DSTL_GA_H_ */
