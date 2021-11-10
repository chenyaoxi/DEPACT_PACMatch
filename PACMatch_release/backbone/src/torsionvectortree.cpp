/*
 * torsionvectortree.cpp
 *
 *  Created on: 2016年12月28日
 *      Author: hyliu
 */

#include "backbone/torsionvectortree.h"

using namespace NSPproteinrep;
using namespace domaintree;

void TorsionVectorTree::init(std::shared_ptr<std::vector<Vector> >vectors,double resolution) {
	torsionvectors_=vectors;
	ndim_=(*vectors)[0].size();
	resolution_=resolution;
	Vector rootcenter(ndim_,0.0);
	Vector halfranges(ndim_,180.0);
	tree_=std::shared_ptr<Tree>(new Tree(Domain(ndim_,rootcenter,halfranges)));
	size_=0.0;
}


