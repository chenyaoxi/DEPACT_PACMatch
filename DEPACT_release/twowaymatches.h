/*
 * twowaymatches.h
 *
 *  Created on: 2018年8月7日
 *      Author: hyliu
 */

#ifndef TWOWAYMATCHES_H_
#define TWOWAYMATCHES_H_
#include<vector>
#include <map>
#include <cassert>
namespace myobcode {
template<typename T1, typename T2>
struct TwoWayMatches{
	TwoWayMatches(){;}
	TwoWayMatches(const T1 *o1, const T2 *o2, const std::vector<int> & matches_in_o2):
	obj1(o1),obj2(o2),seq2(matches_in_o2){
		for(int i=0;i<matches_in_o2.size();++i) {
			seq1.push_back(i);
			map12[i]=matches_in_o2[i];
			map21[matches_in_o2[i]]=i;
		}
	}
	TwoWayMatches(const T1 *o1, const T2 *o2, const std::vector<int> &matches_in_o1,
			const std::vector<int> & matches_in_o2):
	obj1(o1),obj2(o2),seq1(matches_in_o1),seq2(matches_in_o2){
		assert(matches_in_o1.size()==matches_in_o2.size());
		for(int i=0;i<matches_in_o1.size();++i) {
			map12[matches_in_o1[i]]=matches_in_o2[i];
			map21[matches_in_o2[i]]=matches_in_o1[i];
		}
	}
	int o1matchino2(int id_in_o1) const {
		auto find=map12.find(id_in_o1);
		if(find==map12.end()) return -1;
		return find->second;
	}
	int o2matchino1(int id_in_o2) const {
		auto find=map21.find(id_in_o2);
		if(find==map21.end()) return -1;
		return find->second;
	}
	int ithmatchino1(int i) const {return seq1[i];}
	int ithmatchino2(int i) const {return seq2[i];}
	int nmatches() const {
		return map12.size();
	}
	const T1 *obj1 {nullptr};
	const T2 *obj2 {nullptr};
	std::vector<int> seq1;
	std::vector<int> seq2;
	std::map<int,int> map12;
	std::map<int,int> map21;
};
}


#endif /* TWOWAYMATCHES_H_ */
