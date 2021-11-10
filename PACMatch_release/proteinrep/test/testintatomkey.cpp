/*
 * testintatomkey.cpp
 *
 *  Created on: 2016年11月16日
 *      Author: hyliu
 */
#include "proteinrep/intatomkey.h"
#include "proteinrep/pdbrecord.h"
#include "geometry/xyz.h"
#include <iostream>
#include <vector>
#include <algorithm>
using namespace NSPproteinrep;
int main() {
	typedef typename AtomKeyTypeR::Key AtomKey;
	//AtomKey key=AtomKeyTypeR::genKey(0U,"CA",'G',5u,1U);
	std::vector<AtomKey> keys;
	keys.push_back(AtomKeyTypeR::genKey(0U,"CA",0u,"ASP"));
	keys.push_back(AtomKeyTypeR::genKey(8U,"N",1u,"ASPH",2u));
	keys.push_back(AtomKeyTypeR::genKey(9U,"C",1u,"VAL",2u));
	keys.push_back(AtomKeyTypeR::genKey(1U,"O",1u,"PHE",2u));
	std::sort(keys.begin(),keys.end());
	int a=0;
	for( auto key:keys) {
/*		std::cout <<"IntKey: " << key <<std::endl;
	std::cout <<"chain number: " <<AtomKeyTypeR::chainNumber(key) <<std::endl;
	std::cout <<"Position: " << AtomKeyTypeR::posiNumber(key) <<std::endl;
	std::cout <<"Atom Name: " << AtomKeyTypeR::atomName(key) <<std::endl;
	std::cout <<"Residue name " <<AtomKeyTypeR::residueName(key) <<std::endl;
	std::cout <<"Rotamer Code Number: " << AtomKeyTypeR::rotCodeNumber(key) <<std::endl;
	std::cout <<"Is Mainchain? " <<AtomKeyTypeR::isMainChain(key) <<std::endl;
	std::cout <<"Is charged? " <<AtomKeyTypeR::isCharged(key) <<std::endl;
	std::cout <<"Is polar? " <<AtomKeyTypeR::isPolar(key) <<std::endl;
	std::cout <<"Is HB donor? " <<AtomKeyTypeR::isHBDonor(key) <<std::endl;
	std::cout <<"Is HB acceptor? " <<AtomKeyTypeR::isHBAcceptor(key) <<std::endl;
	std::cout <<"Is aromatic? " <<AtomKeyTypeR::isAromatic(key) <<std::endl;
	*/
		PdbRecord record=make_pdbrecord<AtomKeyTypeR>(key,NSPgeometry::XYZ(1.0,2.0,3.0),a++,1);
		std::cout <<record.toString() <<std::endl;
	}
	std::string s{"ATOM      3  C   GLY A 102      64.542  58.454 128.871  1.00 27.48           C"};
	std::cout <<"Line in:" <<std::endl;
	std::cout <<s <<std::endl;
	std::cout <<"Line out:" <<std::endl;
	PdbRecord r(s);
	std::cout <<r.toString()<<std::endl;
}
