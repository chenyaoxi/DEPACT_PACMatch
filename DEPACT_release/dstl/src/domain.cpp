/*
 * domaintree.cpp
 *
 *  Created on: 2016年2月26日
 *      Author: hyliu
 */
#include "dstl/domaintree.h"
using namespace domaintree;
Domain::Domain(int nd, const point_type & cent, const point_type & halves) :
			ndims(nd),center(cent), halfsize(halves) {
		assert(center.size() == ndims);
		assert(halfsize.size() == ndims);
	}
Domain Domain::subdomain(const std::bitset<DOMAINMAXDIM> & code) const {
		point_type newcenter = center;
		point_type newhalves = halfsize;
		for (int d = 0; d < ndims; ++d) {
			if (!code.test(2*d))
				continue;
			newhalves[d] *= 0.5;
			if (code.test(2*d+1))
				newcenter[d] += newhalves[d];
			else
				newcenter[d] -= newhalves[d];
		}
		return Domain(ndims,newcenter, newhalves);
	}
void Domain::print() const {
	std::cout << "Domain:\t";
	for (int d = 0; d < ndims; ++d) {
		std::cout << center[d] << "+-" << halfsize[d] << "\t";
	}
	std::cout << std::endl;
}

