/*
 * sortindex.h
 *
 *  Created on: 2017年3月30日
 *      Author: hyliu
 */

#ifndef DSTL_SORTINDEX_H_
#define DSTL_SORTINDEX_H_
#include <algorithm>
#include <vector>
namespace NSPdstl{
template<typename T>
std::vector<int> sortindex(std::vector<T> & inputvec) {
	std::vector<int> idx;
	for(int i=0;i<inputvec.size();++i) idx.push_back(i);
	auto f=[&inputvec](const int & x,const int & y)->bool{return inputvec[x] < inputvec[y];};
	std::sort(idx.begin(),idx.end(),f);
	return idx;
}

template<typename T1,typename T2>
void sort12(std::vector<T1> &v1, std::vector<T2> &v2) {
	std::vector<int> sortedidx=sortindex(v1);
	std::vector<T1> v1p=v1;
	std::vector<T2> v2p=v2;
	for(int i=0;i<v1p.size();++i){
		v1[i]=v1p[sortedidx[i]];
		v2[i]=v2p[sortedidx[i]];
	}
}

}


#endif /* DSTL_SORTINDEX_H_ */
