/*
 * testcoresegments.cpp
 *
 *  Created on: 2017年5月11日
 *      Author: hyliu
 */


#include "backbone/coresegments.h"
#include "dstl/permutation.h"
#include <set>

using namespace NSPproteinrep;

int main(int argc, char **argv) {
	std::string filename(argv[1]);
	std::ifstream ifs;
	ifs.open(filename.c_str());
	CoreSegments core;
	core.read(ifs);
	for(int s1=0; s1<core.size();++s1) {
		for(int s2=0; s2<core.size(); ++s2) {
			if (s1==s2) continue;
			core.buildclosedloops(s1,s2,1,5);
		}
	}
	for(int s1=0; s1<core.size();++s1) {
		for(int s2=0; s2<core.size(); ++s2) {
			if (s1==s2) continue;
			std::cout <<"Number of loops between segments " << s1 <<
					" and " << s2 << core.loopsbetween(s1,s2)->size() <<std::endl;
		}
	}
	std::vector<int> origin_order;
	for (unsigned int i = 0; i < core.size(); ++i) {
		origin_order.push_back(i);
	}
	int np = NSPdstl::Permutation<int>::factorial(origin_order.size());
	for (unsigned int k = 0; k < np; ++k) {
		std::vector<int> order = NSPdstl::Permutation<int>::getPermutation(
				origin_order, k);
		std::cout << "permutation " << k << ": ";
		std::cout << core.buildlinkedchains(order) <<" possible linkedhains" <<std::endl;
	}
	core.sortlinkedchains();
	int id=0;
	std::set<std::string> writtenchains;
	for(const auto & lc:core.linkedchains()){
		MainChain mc;
		core.buildlinkedmainchain(lc,&mc);
		std::string lcname=core.linkedchainname(lc);
		if(writtenchains.find(lcname) != writtenchains.end()) continue;
		writtenchains.insert(lcname);
		std::string chainname="C"+lcname;
		std::cout <<chainname<< " energy_per_residue: "<<
				core.linkedchainaverageene(lc) <<std::endl;
		std::string pdbname=chainname+".pdb";
		std::string datname=chainname+".dat";
		std::ofstream ofs;
		ofs.open(pdbname.c_str());
		writeSitesToPDB(ofs,mc);
		ofs.close();
		ofs.open(datname.c_str());
		mc.write(ofs);
		ofs.close();
	}
}

