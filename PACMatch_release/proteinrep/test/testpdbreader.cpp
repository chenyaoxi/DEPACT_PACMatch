/*
 * testpdbreader.cpp
 *
 *  Created on: 2016年11月16日
 *      Author: hyliu
 */

#include "proteinrep/pdbreader.h"
#include "proteinrep/aaconformer.h"
#include <iostream>
using namespace NSPproteinrep;

int main() {
	PDBModelID mid;
	mid.pdbid="1eta";
	mid.biounit=2;
	mid.model=2;
	std::string modelpdbfile=extractpdbmodel(mid);
	AAConformersInModel modelconformers;
//	PdbReader reader;
//	reader.readpdb("test.pdb");
	modelconformers.pdbmodelid=mid;
	modelconformers.readpdbfile(modelpdbfile);
//	reader.readpdb(modelpdbfile);
//	std::string chainids=reader.chainids();
	for (int chainnumber=0; chainnumber<modelconformers.conformers.size();++chainnumber){
		std::vector<AAConformer> & conformersinchain=modelconformers.conformers[chainnumber];
		int atomid0=1;
		int resposi=0;
		char chainid=modelconformers.mappdbkeyint->pdbChainID(chainnumber);
		for(auto & c:conformersinchain) {
			auto reskey=modelconformers.mappdbkeyint->pdbResKey(resposi++,chainnumber);
			c.calclocalcrd();
			c.calcglobalcrd();
			std::vector<PdbRecord> pdbrecords=c.make_pdbrecords(chainid,reskey.first,
					reskey.second,atomid0);
			atomid0 += pdbrecords.size();
			for(auto & atm:pdbrecords)
				std::cout << atm.toString()<<std::endl;
		}

	}
	/*
	for(int i=0; i<chainids.size();++i) {
		std::cout << "Chain " << chainids[i] <<" sequence: " <<std::endl;
		const std::vector<std::string> &seq=reader.getaminoacidsequence(chainids[i]);
		for(auto &aatype:seq) std::cout <<" "<< aatype;
		std::cout <<std::endl;
	}
	*/
/*	MapPdbKeyInt mappdbkeyint(reader.records());
	typename PdbReader::RecordsType & records= reader.records();
	for(auto & c: records) {
		for (auto & r:c.second) {
			for (auto & a:r.second) {
				std::cout <<mappdbkeyint.chainNumber(a.chainid) <<"\t"
						<<mappdbkeyint.posiNumber(std::make_pair(a.residueid,a.insertionid),a.chainid)
						<<"\t" <<a.toString() <<std::endl;
			}
		}
	}
	*/
}


