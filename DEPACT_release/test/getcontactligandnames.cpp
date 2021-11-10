/*
 * getcontactligandnames.cpp
 *
 *  Created on: 2018年9月15日
 *      Author: hyliu
 */
#include "analyzecontact.h"
#include "atomcontacts.h"
#include "atomtypessmarts.h"
#include "scorecontact.h"

#include <iostream>
#include <fstream>
using namespace subsitedesign;

int main(int argc,char **argv){
	std::ifstream ifs((std::string(argv[1])));
	while (ifs.good()) {
		int nlatoms;
		std::vector<std::string> atypes;
		std::set<int> excluded;
		std::string molname;
		std::vector<ContactDetails> dts = getdtlsnextmol(ifs, molname,nlatoms, atypes,
				excluded);
		if (dts.empty() || nlatoms == 0)
			continue;
		std::cout <<molname <<std::endl;
	}
}

