/*
 * multitorsiontree.h
 *
 *  Created on: 2016年12月28日
 *      Author: hyliu
 */

#ifndef BACKBONE_TORSIONVECTORTREE_H_
#define BACKBONE_TORSIONVECTORTREE_H_
#include "dstl/domaintree.h"
#include "dstl/domainleaf.h"
#include <memory>
namespace NSPproteinrep{

class TorsionVectorTree {
public:
	typedef domaintree::DomainTree<domaintree::Domain,domaintree::LeafStruct<>> Tree;
	typedef std::vector<double> Vector;
	void init(std::shared_ptr<std::vector<Vector>> vectors,double resolution);
	void insertpoint(long idx){
		std::vector<double> resl(ndim_,resolution_);
		domaintree::LeafOperation<long> leafopt(idx);
		tree_->insertpoint((*torsionvectors_)[idx],resl,leafopt);
		size_+=1.0;
	}
	Tree & gettree(){return *tree_;}
	std::vector<Vector> & gettorsionvectors(){return *torsionvectors_;}
	const std::vector<Vector> & gettorsionvectors() const {return *torsionvectors_;}
	const Tree &gettree() const {return *tree_;}
	double size() {return size_;};
private:
	std::shared_ptr<Tree> tree_;
	std::shared_ptr<std::vector<Vector>> torsionvectors_;
	int ndim_{0};
	double resolution_{0.0};
	double size_;
};

}





#endif /* BACKBONE_TORSIONVECTORTREE_H_ */
