/*
 * testchainpackminimizer.cpp
 *
 *  Created on: 2017年8月17日
 *      Author: hyliu
 */


#include "backbone/chainpackminimizer.h"


using namespace NSPproteinrep;
using namespace NSPdstl;
std::vector<std::vector<BackBoneSite>> extractss(const std::vector<BackBoneSite> &sites){
	std::vector<std::vector<BackBoneSite>> ss;
	bool inss=false;
	char sscode_old='C';
	for(auto &s:sites) {
		char sscode=s.sscode;
		if(sscode == 'H' || sscode =='E'){
			if(inss && sscode ==sscode_old){
				ss.back().push_back(s);
			} else {
				ss.push_back(std::vector<BackBoneSite>());
				ss.back().push_back(s);
				inss=true;
			}
		} else{
			inss=false;
		}
		sscode_old=sscode;
	}
	return ss;
}

int main(int argc, char**argv) {
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
	std::vector<std::vector<BackBoneSite>> ss=extractss(sites);
	std::cout <<"Number of SS elements:" <<ss.size()<<std::endl;
//	writeSitesToPDB(ofs,chain);
	ChainPackMinimizer min;
	MinimizerControl mct;
	mct.maxsteps=500000;
	mct.maxnochangesteps=10000;
	mct.settriphase(0.8,0.4,std::vector<int>({1001,1501,5001}));
	mct.nprint=100;
	min.run(mct,ss);
}
