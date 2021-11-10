/*
 * completechains.cpp
 *
 *  Created on: 2016年12月5日
 *      Author: hyliu
 */

#include <backbone/loopinteractions.h>
#include "dstl/perwalker.h"
#include "dstl/permutation.h"
#include "backbone/linkerselector.h"

#include "backbone/segments.h"
#include <iostream>
#include <sstream>
#include <string>
using namespace NSPproteinrep;

class AtLeaf {
public:
	AtLeaf(LinkerSelector<LinkerInteractions> *s,std::ostream *os): selector_(s), os_(os), chainnumber_(0) {;}
	void operator()(const NSPdstl::PerWalker<std::vector<BackBoneSite> *> & walker) const {
		const std::vector<std::vector<BackBoneSite> *> & selectedlinkers=walker.beads();
		std::vector<BackBoneSite> collection;
		for(unsigned int m=0; m<selector_->segmentseq().size(); ++m) {
			for (auto & s:selector_->getSSsegment(m)) collection.push_back(s);
			if(m+1 <selector_->segmentseq().size())
				for(auto & s: *(selectedlinkers[m])) collection.push_back(s);
		}
		++chainnumber_;
		(*os_) << "Completed Chain " << chainnumber_ << std::endl;
		writeSitesToPDB(*os_,collection);
	}
private:
	LinkerSelector<LinkerInteractions> *selector_;
	std::ostream *os_;
	mutable unsigned int chainnumber_;
};
int main(int argc, char **argv) {
	std::ifstream ifs;
	ifs.open(argv[1]); // segments files
	std::vector<std::vector<BackBoneSite>> segments;
	readsegments(ifs,segments);
	ifs.close();

	LinkerSelector<LinkerInteractions> selector(&segments);
	ifs.open(argv[2]); // linkers
	std::vector<std::vector<BackBoneSite>> linkers;
	while (ifs.good()) {
		char buffer[120];
		ifs.getline(buffer,120);
		std::string str(buffer);
		std::istringstream iss(str);
		unsigned int seg1;
		unsigned int seg2;
		unsigned int length;
		unsigned int ncopies;
		iss>>seg1;
		iss>>seg2;
		iss >>length;
		iss >>ncopies;
		if(!ifs.good()) break;
		for(unsigned int n=0; n<ncopies; ++n) {
			linkers.push_back(std::vector<BackBoneSite>());
			if( !readbackbonesites(ifs,length,linkers.back())){
				std::cout <<"read linkers unsuccessful!" <<std::endl;
				exit(1);
			}
			selector.addLinker(seg1,seg2,&(linkers.back()));
		}
	}
	NSPdstl::PerWalker<std::vector<BackBoneSite> *> walker(segments.size()-1);

	std::ofstream ofs;
	ofs.open("completedchains.pdb");
	unsigned int npermut=NSPdstl::Permutation<unsigned int>::factorial(segments.size());
	AtLeaf atleaf(&selector, &ofs);
	std::vector<unsigned int> original;
	for(unsigned int i=0; i<segments.size(); ++i) original.push_back(i);
	for(unsigned int k=0; k<npermut; ++k) {
		selector.segmentseq()= NSPdstl::Permutation<unsigned int>::getPermutation(original,k);
		walker.enumstep(selector,atleaf);
	}
	ofs.close();
}
