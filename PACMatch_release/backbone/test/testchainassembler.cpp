/*
 * testchainassembler.cpp
 *
 *  Created on: 2017年1月5日
 *      Author: hyliu
 */


#include "backbone/chainassembler.h"
using namespace NSPproteinrep;

int main(int argc, char **argv){
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
	std::ifstream ifs;
	ifs.open(argv[2]); // sites files
	std::vector<std::pair<unsigned int,unsigned int>> loops;
	while (ifs.good()) {
		unsigned int p,l;
		ifs >> p>>l;
		if (ifs.good()) {
			assert(p+l < sites.size());
			loops.push_back(std::make_pair(p,l));
		}
	}
	std::vector<std::vector<BackBoneSite>> elements;
	elements.resize(loops.size()+1, std::vector<BackBoneSite>());
	auto begin=sites.begin();
	unsigned int m=0;
	for(auto & lp:loops) {
		auto end= sites.begin()+lp.first;
		elements[m].resize(end-begin,BackBoneSite());
		std::copy(begin,end, elements[m].begin());
		++m;
		begin=sites.begin()+lp.first+lp.second;
	}
	auto end=sites.end();
	elements[m].resize(end-begin,BackBoneSite());
	std::copy(begin,end, elements[m].begin());
	std::vector<unsigned int> order;
	for(unsigned int i=0;i <elements.size();++i) order.push_back(i);
	std::map<std::pair<unsigned int,unsigned int>,unsigned int> looplengths;
	unsigned int l=0;
	for(auto &lp:loops) {
		looplengths.insert(std::make_pair(std::make_pair(l,l+1),lp.second));
		++l;
	}
	assembleSSelements(elements,
			order,
			looplengths);
}
