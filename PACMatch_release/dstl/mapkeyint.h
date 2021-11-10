/*
 * mapkeyint.h
 *
 *  Created on: 2016年11月17日
 *      Author: hyliu
 */

#ifndef DSTL_MAPKEYINT_H_
#define DSTL_MAPKEYINT_H_
#include <map>
#include <vector>
namespace NSPdstl {

template <typename KEY>
class MapKeyInt {
public:
	typedef KEY KeyType;
	MapKeyInt(){;}
	template <typename VAL>
	void init (const std::map<KeyType,VAL> & map) {
		for(auto & e:map){
			if(keynumbers_.find(e.first) == keynumbers_.end()){
				keys_.push_back(e.first);
				keynumbers_.insert(std::make_pair(e.first, keys_.size()-1));
			}
		}
	}
	KeyType key(unsigned int i) const {return keys_.at(i);}
	unsigned int keynumber (const KeyType & key) const {
		if (keynumbers_.count(key) == 0)
			return -10;
		else
			return keynumbers_.at(key);
	}
	template <typename VAL>
	VAL mapvalue(const std::map<KeyType,VAL> & map,unsigned int i) const {
		return map.at(keys_(i));
	}
private:
	std::vector<KeyType> keys_;
	std::map<KeyType,unsigned int> keynumbers_;
};
}



#endif /* DSTL_MAPKEYINT_H_ */
