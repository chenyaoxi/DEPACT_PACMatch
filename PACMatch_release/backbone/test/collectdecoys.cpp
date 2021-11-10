/*
 * evaluatedecoys.cpp
 *
 *  Created on: 2017年8月9日
 *      Author: hyliu
 */

#include "backbone/backboneenergy.h"
#include "backbone/backbonemoves.h"
#include "dstl/randomengine.h"
#include <iostream>
#include <fstream>
#include <sstream>
using namespace NSPproteinrep;
void writedecoy(std::ostream &os, int decoyid,
		const std::vector<BackBoneSite> &chain) {
	os << decoyid << " " << chain[0].pdbid << " " << chain.size() << std::endl;
	for (auto s : chain) {
		if (s.sscode == ' ')
			s.sscode = 'x';
		if (s.resname == "CISPRO")
			s.resname = "PRO";
		os << s.toString();
	}
}
bool readnextdecoy(std::istream &is, int *decoyid,
		std::vector<BackBoneSite> *chain) {
	std::string pdbid;
	int chainsize;
	char buffer[120];
	int nsection = 0;
	is.getline(buffer, 120);
	if (!is.good())
		return false;
	std::string line = std::string(buffer);
	std::stringstream sstr(line);
	sstr >> *decoyid >> pdbid >> chainsize;
//	if(!sstr.eof() && !sstr.good()) return false;
	chain->clear();
	return readbackbonesites(is, chainsize, *chain);
}
std::set<int> sspositions(const std::vector<BackBoneSite> &chain) {
	std::set<int> res;
	for (int i = 0; i < chain.size(); ++i) {
		char ss = chain.at(i).sscodechar();
		if (ss == 'H' || ss == 'E')
			res.insert(i);
	}
	return res;
}
int main(int argc,char **argv){
	std::ifstream ifs;
	ifs.open("alldecoys.dat");
	initenergycontrols();
	BackBoneEnergy energy;
	std::ofstream ofs;
	ofs.open("decoyenergies.dat");
	while (true){
		std::vector<BackBoneSite>chain;
		int decoyid;
		if(!readnextdecoy(ifs, &decoyid,&chain)) break;
		BackBoneEnergy::preparechain(&chain);
		EnergyComponents ecomp;
		energy.totalenergy(chain,&ecomp);
		ofs <<chain[0].pdbid<<" "<< decoyid;
		for(auto e:ecomp.energies) ofs <<" "<<e;
		ofs <<std::endl;
	}
}
/*
int main(int argc, char**argv) {
	std::vector<std::string> directories{"random2347", "random345",   "random55",
		"random773",  "random851",
		"random131",   "random2467",  "random3571",  "random7361",  "random8331",
		"random9231","random2335","random33","random3759", "random7717", "random8453",
		"random9453"};
	std::vector<std::ifstream> ifss;
	for(auto &dir:directories){
		std::string filename=dir+"/"+"decoys.dat";
		ifss.push_back(std::ifstream());
		ifss.back().open(filename.c_str());
	}
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
	std::ofstream ofs;
	ofs.open("alldecoys.dat");
	int maxchains = 100;
	int nchains = 0;
	auto itstart = sites.begin();
	int ndecoys = 0;
	int decoysperchain = 5;
	while (nchains < maxchains && itstart < sites.end()) {
		auto itend = itstart;
		while (!chainendsite(itend, sites.end()))
			++itend;
		int chainlength = itend - itstart + 1;
		if (chainlength >= 150 && chainlength < 500) {
			if (fragstartsite(itstart, sites.end(), chainlength, std::string(),
					false)) {
				std::vector<BackBoneSite> chain(chainlength);
				std::copy(itstart, itstart + chainlength, chain.begin());
//				std::string pdbfile = chain[0].pdbid + "_ref" + ".pdb";
//				std::ofstream ofspdb;
//				ofspdb.open(pdbfile.c_str());
//				writeSitesToPDB(ofspdb, chain);
//				ofspdb.close();
				writedecoy(ofs, ndecoys++,chain);
				for(auto &is:ifss){
					for(int i=0;i<5;++i){
						std::vector<BackBoneSite> decoy;
						int decoyid;
						if(readnextdecoy(is,&decoyid,&decoy)){
							writedecoy(ofs,ndecoys++,decoy);
						} else {
							break;
						}
					}
				}
				++nchains;
			}
		}
		itstart = itend + 1;
	}
	ofs.close();
}
*/


