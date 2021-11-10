/*
 * smartsfragments.h
 *
 *  Created on: 2018年8月6日
 *      Author: hyliu
 */

#ifndef SMARTSFRAGMENTS_H_
#define SMARTSFRAGMENTS_H_
#include "twowaymatches.h"
#include <map>
#include <string>
#include <vector>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <memory>
namespace myobcode{

struct SMARTSFragment{
	std::string smarts;   //smarts string
	std::vector<int> scatoms; //atoms for specific contacts
	std::vector<std::string> scatomtypes;
	static const std::map<std::string,SMARTSFragment> & getfragments
		(const std::string &filename="smartsfragments.txt");
};

}

#endif /* SMARTSFRAGMENTS_H_ */
