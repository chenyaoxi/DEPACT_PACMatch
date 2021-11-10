/*
 * bfdecoupling.cpp
 *
 *  Created on: 2019年1月24日
 *      Author: yxchen
 */

#include "bfdecoupling.h"
using namespace OpenBabel;
using namespace subsitedesign;

std::string myobcode::bfdecoupling(
		const BasicFragment *sf, OpenBabel::OBMol *mol) {
	std::string out = "Nan";
	OpenBabel::OBSmartsPattern smarts;
	smarts.Init(sf->smarts);
	std::vector<std::vector<int> > maplist;
	if (smarts.Match(*mol)) {
		maplist = smarts.GetMapList();
		if (maplist.size() > 0) out = "Yes";
	}
	return out;
}
