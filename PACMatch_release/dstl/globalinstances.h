/*
 * globalinstances.h
 *
 *  Created on: 2017年8月19日
 *      Author: hyliu
 */

#ifndef DSTL_GLOBALINSTANCES_H_
#define DSTL_GLOBALINSTANCES_H_

#include <map>
#include <string>
namespace NSPdstl {
template <typename OBJ>
class GlobalInstances{
public:
	static OBJ & getstaticinstance(){
		static OBJ obj;
		return obj;
	}
	static OBJ &getnamedinstance(const std::string & id){
		static std::map<std::string, OBJ>  namedinstances;
		if(id.empty()) return getstaticinstance();
		auto it=namedinstances.find(id);
		if(it == namedinstances.end()){
#ifdef _OPENMP
#pragma omp critical(namedinstance_global)
			{
				auto itx=namedinstances.find(id);
				if(it ==namedinstances.end()){
#endif
			namedinstances.insert(std::make_pair(id,OBJ()));
#ifdef _OPENMP
#pragma omp flush
				}
			}
#endif
			return namedinstances.at(id);
		} else {
			return it->second;
		}
	}
	static OBJ &copystaticinstance(const std::string &id){
		getnamedinstance(id)=getstaticinstance();
		return getnamedinstance(id);
	}
	static OBJ &copynamedinstance(const std::string &from, const std::string &to){
		getnamedinstance(to)=getnamedinstance(from);
		return getnamedinstance(to);
	}
	static OBJ &getinstance(const std::string &id=std::string()){
		if(id.empty()) return getstaticinstance();
		else return getnamedinstance(id);
	}
};
}

#endif /* DSTL_GLOBALINSTANCES_H_ */
