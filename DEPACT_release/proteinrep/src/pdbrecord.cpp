/*
 * pdbatomrecord.cpp
 *
 *  Created on: 2016年11月16日
 *      Author: hyliu
 */
#include "proteinrep/pdbrecord.h"
#include "dataio/inputlines.h"
#include <boost/algorithm/string.hpp>
using namespace NSPproteinrep;

PdbRecord::PdbRecord(const std::string &line) {
	const std::vector<int> fw { 6, 6, 2, 2, 1, 4, 1, 4, 1, 11, 8, 8, 6, 6, 10, 2 };
	std::string eline =line+std::string(30,' ');
	init(NSPdataio::parseline(eline, fw));
}

void PdbRecord::init(const std::vector<std::string> &line) {
	label = boost::trim_copy(line[LABEL]);
	namesymbol = boost::trim_copy(line[NAMESYMBOL]);
	namemodifier = boost::trim_copy(line[NAMEMODIFIER]);
	atomname = namesymbol + namemodifier;
	if(!line[CONFORMERID].empty())
		conformerid=line[CONFORMERID][0];
	if(!line[INSERTIONID].empty()){
		insertionid=line[INSERTIONID][0];
	}
	chainid = line[CHAINID][0];
	if (chainid == ' ')
		chainid = CHAINID_DEFAULT;
	residuename = line[RESIDUENAME];
	try {
		residueid = std::stoi(line[RESIDUEID]);
		atomid = std::stol(line[ATOMID]);
		x = std::stod(line[X]);
		y = std::stod(line[Y]);
		z = std::stod(line[Z]);
	} catch (std::exception &e) {
		std::string mesg =
				"Error converting to numeric when initialize a PdbAtomRecord\n";
		std::cout << mesg << std::endl;
		failure = true;
//		NSPdataio::attendError(mesg);
	}
	try {
		occupation = std::stod(line[OCCUPATION]);
		bfactor = std::stod(line[BFACTOR]);
		if(line[ELEMENTNAME].size() <2 ) elementname[1]=line[ELEMENTNAME][0];
		else {
			elementname[0] = line[ELEMENTNAME][0];
			elementname[1] = line[ELEMENTNAME][1];
		}
	} catch (std::exception &e) {
		;
	}
}

std::string PdbRecord::toString() const {
	// this format start position of aotmname not standard for 2-letter chemical symbols
	const std::string FMT {
			"%-6s%5d %2s%-2s%c%3s %1c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %c%c" };
	char s[100];
	if (sprintf(s, FMT.c_str(), label.c_str(), atomid, namesymbol.c_str(),
			namemodifier.c_str(), conformerid, residuename.substr(0,3).c_str(), chainid,
			residueid, insertionid, x, y, z, occupation, bfactor, elementname[0],elementname[1])
			< 0)
		NSPdataio::attendError("Error printing PdbRecord");
	return (std::string(s));
}

