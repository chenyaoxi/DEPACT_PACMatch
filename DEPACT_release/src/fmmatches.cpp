/*
 * fmmatches.cpp
 *
 *  Created on: 2018年8月11日
 *      Author: hyliu
 */
#include "fmmatches.h"
using namespace OpenBabel;
using namespace subsitedesign;
using namespace myobcode;
std::vector<std::shared_ptr<FMMatches>> myobcode::findfmmatches(
		const BasicFragment *sf, OpenBabel::OBMol *mol, bool unique) {
	OpenBabel::OBSmartsPattern smarts;
	smarts.Init(sf->smarts);
	std::vector<std::vector<int> > maplist;
	std::vector < std::shared_ptr < FMMatches >> res;
	if (smarts.Match(*mol)) {
		if(unique) maplist = smarts.GetUMapList();
		else maplist=smarts.GetMapList();
		for (auto itr = maplist.begin(); itr != maplist.end(); itr++) {
			res.push_back(
					std::shared_ptr < FMMatches
							> (new FMMatches(sf, mol, *itr)));
		}
	}
	return res;
}
/*std::vector<std::string> myobcode::findatomtypes(const std::map<std::string,std::vector<std::shared_ptr<FMMatches>>> &fmmatches){
	const OBMol * mol;
	for(auto & fv:fmmatches){
		if(fv.second.empty()) continue;
		mol=fv.second.at(0)->obj2;
		break;
	}
	int natoms=mol->NumAtoms();
	std::vector<std::string> atomtypes(natoms,"unspecified");
	for(auto & fv:fmmatches){
		if(fv.second.empty()) continue;
		const BasicFragment* bfg=fv.second.at(0)->obj1;
		const std::vector<int> & spatoms=bfg->spatoms;
		const std::vector<std::string> & sptypes=bfg->sptypes;
		for(auto &fm:fv.second){
			const std::vector<int> & seq2=fm->seq2;
			for(int i=0; i<spatoms.size();++i){
				int j=fm->o1matchino2(spatoms.at(i))-1; //atoms in OBMol starts from 1
				atomtypes[j]=sptypes[i];
			}
		}
	}
	return atomtypes;
}*/



