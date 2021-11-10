/*
 * readpdbchains.cpp
 *
 *  Created on: 2017年6月27日
 *      Author: hyliu
 */

#include "dataio/pdbfolder.h"
#include "dataio/inputlines.h"
#include "proteinrep/aaconformer.h"
#include <string>
using namespace NSPdataio;
using namespace NSPproteinrep;

int main(int argc, char **argv){
	TextLines txts;
	txts.init(std::string(argv[1]));
	std::multimap<std::string,char> pdbchainids;
	for(auto & word:txts.lines()){
		pdbchainids.insert(
				std::make_pair(word.substr(0,4),word[4]));
	}
	std::string pdbroot="";
	if(argc>2)pdbroot=std::string(argv[2]);
	PdbFolder &pdbfolder= PdbFolder::getinstance(pdbroot);
	std::vector<std::vector<AAConformer>> allconformers;
	for(auto & pdbchain:pdbchainids){
		std::string pdbid=pdbchain.first;
		char chainid=pdbchain.second;
		std::cout <<"Reading " <<pdbid << " chain " <<chainid<<std::endl;
		std::string localfilename=pdbfolder.getfile(pdbid);
		if(localfilename.empty()) {
			std::cout <<"Get PDB file for " << pdbid <<"failed." <<std::endl;
			abort();
		}
		PDBModelID mid;
		mid.pdbid=pdbid;
		AAConformersInModel modelconformers;
		modelconformers.pdbmodelid=mid;
		modelconformers.readpdbfile(localfilename);
		if (chainid=='A') chainid=' ';
		const std::vector<AAConformer> & conformers=modelconformers.conformersinchain(chainid);
		allconformers.push_back(std::vector<AAConformer>());
		std::vector<AAConformer> & conf=allconformers.back();
		conf.resize(conformers.size());
		std::copy(conformers.begin(),conformers.end(),conf.begin());
	}
}


