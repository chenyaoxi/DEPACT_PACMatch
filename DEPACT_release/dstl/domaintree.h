/*
 * domaintree.h
 *
 *  Created on: 2016年2月26日
 *      Author: hyliu
 */

#ifndef DOMAINTREE_H_
#define DOMAINTREE_H_

#include "dstl/nnearest.h"
#include <vector>
#include <map>
#include <cassert>
#include <bitset>
#include <iostream>
#include <algorithm>
#include <memory>

#define DOMAINMAXDIM 64

namespace domaintree {

struct Domain {
	typedef std::vector<double> point_type;
	int ndims;
	point_type center;
	point_type halfsize;
	Domain():ndims(0) {;}
	Domain(int nd, const point_type & cent, const point_type & halves);
	Domain subdomain(const std::bitset<DOMAINMAXDIM> & code) const;
	template<typename T>
	bool contain(const T & point) const {
//		assert(point.size() == ndims);
		for (int d = 0; d < ndims; ++d) {
			if (point[d] < center[d] - halfsize[d]-0.00000001
					|| point[d] > center[d] + halfsize[d]+0.0000001)
				return false;
		}
		return true;
	}
	void print() const;
};
/*
 template<typename T, typename DOMAIN_>
 double point_domain_dist2(T point, const DOMAIN_ & domain) {
 double dist2;
 for (int d = 0; d < DOMAIN_::ndims; ++d) {
 double tmp = point[d] - domain.center[d];
 dist2 += tmp * tmp;
 }
 return dist2;
 }
 */
template<typename T, typename D, typename D2LEAF>
double point_domain_mindist2(T point, const D& domain, D2LEAF *d2leaf) {
	typedef typename D2LEAF::crd_type crd_type;
//	if (domain.contain(point))
//		return 0.0;
	double dist2{0.0};
	std::vector<double> tmp=crd_type::diff(point,domain.center);
	for (int d = 0; d < domain.ndims; ++d) {
//		double tmp = crd_type::diff(point[d], domain.center[d]);
		if(tmp[d] <0) tmp[d]=-tmp[d];
		if( tmp[d] < domain.halfsize[d])
			tmp[d]=0;
		else
			tmp[d] -= domain.halfsize[d];
		dist2 += tmp[d] * tmp[d];
	}
	return dist2;
}

template<typename DOMAIN_, typename D2LEAF>
double inter_domain_dist2(const DOMAIN_ & d1, const DOMAIN_ & d2,
		D2LEAF *d2leaf) {
	typedef typename D2LEAF::crd_type crd_type;
	double dist2{0.0};
	std::vector<double> tmp =crd_type::diff(d1.center, d2.center);
	for (int d = 0; d < d1.ndims; ++d) {
//		double tmp = crd_type::diff(d1.center[d], d2.center[d]);
		dist2 += tmp[d] * tmp[d];
	}
	return dist2;
}

template<typename DOMAIN_, typename LEAF>
class DomainTree {
public:
	typedef typename DOMAIN_::point_type point_type;
	DomainTree<DOMAIN_,LEAF>(): ndims_(0){;}
	DomainTree<DOMAIN_, LEAF>(const DOMAIN_ &d) :
			leaf_(), domain_(d) {
		ndims_ = d.ndims;
	}
	DOMAIN_ & domain() {
		return domain_;
	}
	const DOMAIN_ & domain() const {
		return domain_;
	}
	template<typename LEAFOP>
	void insertpoint(const point_type & point, const point_type & resolution,
			LEAFOP &leafop) {
		/*		domain_.print();
		 for(auto x:point)
		 std::cout <<"\t" <<x;
		 std::cout <<std::endl;
		 */
		assert(domain_.contain(point));
/*		if( !domain_.contain(point)) {
			std::cout <<"ndims_"<<ndims_;
			for(int d=0; d<ndims_;d++) {
					std::cout
						<<point[d] <<" ";
				}
			std::cout <<std::endl;
			abort();
		}
		*/
		std::vector<bool> torefine(ndims_);
		if (reach_resolution(resolution, torefine)) {
			leafop(domain_, leaf_);
			return;
		}
		std::bitset<DOMAINMAXDIM> code;
		code.reset();
		for (int d = 0; d < ndims_; ++d) {
			if (torefine[d]) {
				code.set(2 * d, 1);
				if (point[d] >= domain_.center[d])
					code.set(2 * d + 1, 1);
			}
		}
		std::string scode=code.to_string();
		auto found = subdomains_.find(scode);
		if (found != subdomains_.end()) {
			found->second->insertpoint(point, resolution, leafop);
		} else {
			DOMAIN_ sub = domain_.subdomain(code);
			DomainTree<DOMAIN_, LEAF> * subtree = new DomainTree<DOMAIN_, LEAF>(
					sub);
			subdomains_.insert(std::make_pair(scode, subtree));
			subtree->insertpoint(point, resolution, leafop);
		}
	}
	bool reach_resolution(const point_type & resolution,
			std::vector<bool> & torefine) {
		bool res = true;
		for (int d = 0; d < ndims_; d++) {
			torefine[d] = false;
			if (domain_.halfsize[d] > resolution[d]) {
				torefine[d] = true;
				res = false;
			}
		}
		return res;
	}

	template<typename D2LEAF>
	void findneighbor(const point_type & point, D2LEAF &d2leaf,
			double &disbound) {
		if (subdomains_.empty()) {
//			d2leaf(leaf_, disbound); //update actual nearest distance to a leaf
			d2leaf(leaf_,point,disbound);
		}
		std::vector<std::pair<std::string, double>> keyd2pairs;
		sortsubdomains(point, keyd2pairs, &d2leaf);
		for (auto & p : keyd2pairs) {
			if (disbound < 0 || p.second < disbound) {
				subdomains_.at(p.first)->findneighbor(point, d2leaf, disbound);
			}
		}
	}
	~DomainTree() {
		for (auto &map : subdomains_)
			delete map.second;

	}
	std::shared_ptr<LEAF> & leaf() {
		return leaf_;
	}
private:
	DOMAIN_ domain_;
	std::map<std::string, DomainTree *> subdomains_;
	std::shared_ptr<LEAF> leaf_;
	int ndims_;
	template<typename D2LEAF>
	void sortsubdomains(const point_type & point,
			std::vector<std::pair<std::string, double>> & keyd2pairs,
			D2LEAF *d2leaf) const {
		keyd2pairs.clear();
		for (auto & elem : subdomains_) {
			double d2=point_domain_mindist2(point, elem.second->domain(),
									d2leaf);
			keyd2pairs.push_back(std::make_pair(elem.first,d2));
		}
		std::sort(keyd2pairs.begin(), keyd2pairs.end(),
				[](std::pair<std::string,double> e1,std::pair<std::string,double> e2)
				->bool {return e1.second < e2.second;});
	}
};
}

#endif /* DOMAINTREE_H_ */
