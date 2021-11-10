/*
 * stlutil.h
 *
 *  Created on: 2017年8月3日
 *      Author: hyliu
 */

#ifndef DSTL_STLUTIL_H_
#define DSTL_STLUTIL_H_
#include <string>
#include <vector>
#include <map>
namespace NSPdstl{
template<typename KEY,typename VAL>
std::vector<KEY> getkeyvec(const std::map<KEY,VAL> &map){
	std::vector<std::string> res;
	for(auto &t:map) res.push_back(t.first);
	return res;
};
}


#endif /* DSTL_STLUTIL_H_ */
