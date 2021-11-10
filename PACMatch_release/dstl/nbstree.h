/*
 * nbstree.h
 *
 *  Created on: 2018年6月27日
 *      Author: hyliu
 */

#ifndef NBSTREE_H_
#define NBSTREE_H_
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <cassert>
#include <iostream>
namespace NSPdstl {
struct NBSTreePar{
	std::vector<double> startvals;
	std::vector<double> binwidths;
	int ndim{0};
};
struct NBSQuery{
	std::vector<double> vals;
	std::vector<double> lcuts;
	std::vector<double> rcuts;
};
template<typename POINT>
class NBSNode{
public:
	NBSNode (int dim=0,NBSTreePar *tp=nullptr):dim_(dim),treepar_(tp){;}
	std::pair<std::set<int>::const_iterator, std::set<int>::const_iterator>
	neighborrange(double centerval,double lcut,double rcut) const{
		int lbinid=binid(centerval-lcut);
		int rbinid=binid(centerval+rcut);
		std::set<int>::const_iterator itl=binset_.lower_bound(lbinid);
		std::set<int>::const_iterator itr=binset_.upper_bound(rbinid);
		return std::make_pair(itl,itr);
	}
	void addpoint(POINT p, std::vector<double> vals){
		if (dim_>=treepar_-> ndim){
			if(!points_) points_=std::shared_ptr<std::vector<POINT>>(new std::vector<POINT>());
			points_->push_back(p);
			++branchsize_;
		} else {
			++branchsize_;
			int bin=binid(vals[dim_]);
			if(binset_.find(bin) == binset_.end()){
				binset_.insert(bin);
				children_[bin]=std::shared_ptr<NBSNode>(new
						NBSNode<POINT>(dim_+1,treepar_));
			}
			children_[bin]->addpoint(p,vals);
		}
	}
	void findneighbors(const NBSQuery &query,
			std::vector<POINT> *neighbors) const {
		if (dim_>=treepar_-> ndim){
			for( auto & p:*points_) neighbors->push_back(p);
		} else {
			auto rng=neighborrange(query.vals[dim_],query.lcuts[dim_],
					query.rcuts[dim_]);
			for(auto it=rng.first; it!=rng.second; ++it){
				children_.at(*it)->findneighbors(query,neighbors);
			}
		}
	}
private:
	std::set<int> binset_;
	std::map<int,std::shared_ptr<NBSNode>> children_;
	NBSTreePar *treepar_;
	int dim_;
	int binid(double val) const {
		if(val<treepar_->startvals[dim_]) return 0;
		return (val-treepar_->startvals[dim_])/treepar_->binwidths[dim_];
	}
	int branchsize_{0};
	std::shared_ptr<std::vector<POINT>> points_;
};
template<typename POINT>
class NBSTree{
public:
	NBSTree(const std::vector<double> &startvals,
			const std::vector<double> &binwidths):
			root(0,&treepar){
		treepar.startvals=startvals;
		treepar.binwidths=binwidths;
		treepar.ndim=startvals.size();
	}
	void addpoint(POINT p, const std::vector<double> &valvec){
		root.addpoint(p,valvec);
	}
	void findneighbors(const NBSQuery & query,
			std::vector<POINT> *neighbors){
		root.findneighbors(query,neighbors);
	}
private:
	NBSTreePar treepar;
	NBSNode<POINT> root;
};
}



#endif /* NBSTREE_H_ */
