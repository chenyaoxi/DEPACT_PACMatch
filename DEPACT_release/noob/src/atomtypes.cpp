/*
 * atomtypes.cpp
 *
 *  Created on: 2018年8月16日
 *      Author: hyliu
 */
#include "atomtypes.h"
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <dataio/datapaths.h>
using namespace NSPdataio;
using namespace subsitedesign;
const std::map<std::string, AtomType> & AtomType::getmap(
		const std::string &filename) {
	static std::map<std::string, AtomType> map;
	static bool initialized { false };
	if (!initialized) {
		initmap(getenvpath("DEPACT_DATAPATH")+filename, map);
		initialized = true;
	}
	return map;
}
std::string subsitedesign::proteinatomtype(const std::string &residue, std::string atom){
	static const std::set<std::string> MCATOMS{"N","CA","C","O"};
	if(atom=="N" && residue =="PRO") return "PRO_N";
	if(atom=="CA" && residue=="GLY") return "GLY_CA";
	if(atom=="OXT") return "OXT";
	if(MCATOMS.find(atom) != MCATOMS.end()) return std::string("MC_")+atom;
	if(residue=="PHE" || residue=="TYR"){
		if(atom=="CD1"||atom=="CD2") atom="CD";
		else if(atom=="CE1"||atom=="CE2") atom="CE";
	}else if(residue=="ASP"){
		if(atom=="OD1"||atom=="OD2") atom="OD";
	} else if(residue=="GLU"){
		if(atom=="OE1" ||atom=="OE2") atom="OE";
	} else if(residue=="ARG"){
		if(atom=="NH1" ||atom=="NH2") atom="NH";
	}
	return residue+"_"+atom;
}
std::string subsitedesign::ptype2(const std::string &ptype1){
	static std::map<std::string,std::string> type2map;
	static bool first{true};
	if(first){
		first=false;
		std::ifstream ifs("proteinatomtypes.txt");
		while(ifs.good()){
			std::string line;
			std::getline(ifs,line);
			if(!ifs.good()) break;
			std::stringstream sstr(line);
			std::string name1;
			std::string name2;
			sstr >>name1 >>name2;
			type2map[name1]=name2;
		}
	}
	if(type2map.find(ptype1)== type2map.end()) return std::string();
	return type2map[ptype1];
}
std::string subsitedesign::proteinatomtype2(const std::string &residue, std::string atom){
	std::string ptype1=proteinatomtype(residue,atom);
	return ptype2(ptype1);
}
void AtomType::initmap(const std::string &filename,
		std::map<std::string, AtomType> &map) {
	std::ifstream ifs(filename);
	typedef boost::tokenizer<boost::char_separator<char> > Tokenizer_char;
	typedef boost::tokenizer<boost::offset_separator> Tokenizer_int;
	if (!ifs.good()) {
		std::cout << "Cannot open atomtypes file " << filename << std::endl;
		exit(1);
	}
	std::vector<std::vector<std::string>> inputlines;
	char inputline[500];
	boost::char_separator<char> sep(" \t");/*! char delimiter set is  " \t"*/
	while (ifs.getline(inputline, 500)) {
		if (inputline[0] == '#')
			continue;
		std::string line(inputline);
		Tokenizer_char tok(line, sep);
		if (tok.begin() == tok.end())
			continue;
		inputlines.push_back(std::vector<std::string>());
		for (auto it = tok.begin(); it != tok.end(); ++it)
			inputlines.back().push_back(std::string(*it));
	}
	int natomtypes;
	for (int lidx = 0; lidx < inputlines.size(); ++lidx) {
		std::vector<std::string> &words = inputlines[lidx];
		if (lidx == 0) {
			natomtypes = std::stoi(words[0]);
		} else if (lidx <= natomtypes) {
			if(words.size() !=(AtomType::NFEATURES + 3)){
				std::cout <<words[2]<<std::endl;
			}
			assert(words.size() == (AtomType::NFEATURES + 3));
			AtomType::AtomFeatures f(0);
			f.set(0); //specified type
			for (int i = 4; i < AtomType::NFEATURES+3; ++i) {
				if (std::stoi(words[i]) != 0)
					f.set(i-3);
			}
			map.insert(
					std::make_pair(words[2], AtomType(words[0], words[2],words[3], f)));
			map[words[2]].codename0=words[1]+"_";
		}
	}
	map["Un0"] = AtomType("Un0","Un0");
	map["Un0"].codename0="Un0";
}
std::string subsitedesign::getcodename0(const std::string &atype){
	return AtomType::getmap().at(atype).codename0;
}
std::string subsitedesign::getcodename(const std::string &atype){
	return AtomType::getmap().at(atype).codename;
}
