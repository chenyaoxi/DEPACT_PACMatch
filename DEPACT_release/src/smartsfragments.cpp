/*
 * smartfragments.cpp
 *
 *  Created on: 2018年8月6日
 *      Author: hyliu
 */

#include "smartsfragments.h"
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <set>
using namespace myobcode;

typedef boost::tokenizer<boost::char_separator<char> > Tokenizer_char;
typedef boost::tokenizer<boost::offset_separator> Tokenizer_int;
const std::map<std::string, SMARTSFragment> & SMARTSFragment::getfragments(
		const std::string &filename) {
	static std::map<std::string, SMARTSFragment> map;
	static bool initialized { false };
	if (!initialized) {
		initialized = true;
		std::ifstream ifs;
		ifs.open(filename.c_str());
		if (!ifs.good()) {
			std::cout << "Cannot open SMARTS_fragments file " << filename
					<< std::endl;
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
		for (int lidx = 0; lidx < inputlines.size(); ++lidx) {
			std::vector<std::string> &words = inputlines[lidx];
			SMARTSFragment sf;
			sf.smarts = words[1];
			for (int i = 2; i < words.size(); ++i)
				sf.scatoms.push_back(std::stoi(words[i]));
			map[words[0]] = sf;
		}
	}
	return map;
}
std::vector<std::shared_ptr<FMMatches>> myobcode::findfmmatches(
		const SMARTSFragment *sf, OpenBabel::OBMol *mol, bool unique) {
	OpenBabel::OBSmartsPattern smarts;
	smarts.Init(sf->smarts);
	std::vector<std::vector<int> > maplist;
	std::vector < std::shared_ptr < FMMatches >> res;
	if (smarts.Match(*mol)) {
		if(unique) maplist = smarts.GetUMapList();
		else maplist=smarts.GetMapList();
		for (auto itr = maplist.begin(); itr != maplist.end(); itr++) {
			res.push_back(
					std::shared_ptr < FMMatches
							> (new FMMatches(sf, mol, *itr)));
		}
	}
	return res;
}

