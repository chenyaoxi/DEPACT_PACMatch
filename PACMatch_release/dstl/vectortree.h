/*
 * vectortree.h
 *
 *  Created on: 2017年12月5日
 *      Author: hyliu
 */

#ifndef DSTL_VECTORTREE_H_
#define DSTL_VECTORTREE_H_
#include "dstl/domaintree.h"
#include "dstl/domainleaf.h"
#include <memory>
namespace domaintree {
class VectorTree {
public:
	typedef DomainTree<Domain,LeafStruct<>> Tree;
	typedef std::vector<double> Vector;
	void init(std::shared_ptr<std::vector<Vector>> vectors,double resolution,double vmin,
			double vmax){
		vectors_=vectors;
		ndim_=(*vectors)[0].size();
		resolution_=resolution;
		double halfrange=0.5*(vmax-vmin);
		double center=vmin +halfrange;
		Vector rootcenter(ndim_,center);
		Vector halfranges(ndim_,halfrange);
		tree_=std::shared_ptr<Tree>(new Tree(Domain(ndim_,rootcenter,halfranges)));
		size_=0.0;
	}
	void insertpoint(long idx){
		std::vector<double> resl(ndim_,resolution_);
		domaintree::LeafOperation<long> leafopt(idx);
		tree_->insertpoint((*vectors_)[idx],resl,leafopt);
		size_+=1.0;
	}
	Tree & gettree(){return *tree_;}
	std::vector<Vector> & getvectors(){return *vectors_;}
	const std::vector<Vector> & getvectors() const {return *vectors_;}
	const Tree &gettree() const {return *tree_;}
	double size() {return size_;};
private:
	std::shared_ptr<Tree> tree_;
	std::shared_ptr<std::vector<Vector>> vectors_;
	int ndim_{0};
	double resolution_{0.0};
	double size_;
};
}


#endif /* DSTL_VECTORTREE_H_ */
