/*
 * basicfragment.cpp
 *
 *  Created on: 2018年8月9日
 *      Author: hyliu
 */

#ifndef BASICFRAGMENT_CPP_
#define BASICFRAGMENT_CPP_
#include "basicfragment.h"

/*
 * smartfragments.cpp
 *
 *  Created on: 2018年8月6日
 *      Author: hyliu
 */

#include "basicfragment.h"
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>
#include <set>
using namespace subsitedesign;

void BasicFragment::initmap(std::map<std::string, BasicFragment> *frgmap,
		const std::string &filename) {
		typedef boost::tokenizer<boost::char_separator<char> > Tokenizer_char;
		typedef boost::tokenizer<boost::offset_separator> Tokenizer_int;
		std::string file = filename; //TODO: add environmental path
		std::ifstream ifs;
		ifs.open(file.c_str());
		if (!ifs.good()) {
			std::cout << "Cannot open basicfragments file " << file
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
				BasicFragment fg;
				fg.smarts = words[1];
				(*frgmap)[words[0]] = fg;
			}
}

#endif /* BASICFRAGMENT_CPP_ */
