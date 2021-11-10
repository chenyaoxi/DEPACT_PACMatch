/*
 * targetbf.cpp
 *
 *  Created on: 2019年1月24日
 *      Author: yxchen
 */
#include "basicfragment.h"
#include <iostream>
#include "bfdecoupling.h"
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
using namespace subsitedesign;
using namespace OpenBabel;
using namespace myobcode;

int main(int argc, char **argv) {
//target.sdf
	const std::map<std::string, BasicFragment> &map =
			BasicFragment::getmap();
	OBConversion obconversion;
	OBMol targetmol;
	obconversion.SetInFormat("sdf");
	std::string targetsdf(argv[1]);
	obconversion.ReadFile(&targetmol, targetsdf);
	std::ofstream ofs;
	std::ifstream ifs;
	ifs.open(targetsdf.c_str());
	if (!ifs.good()) {
		std::cout << "sdf failure" << std::endl;
		exit(1);
	}
	std::string line;
	line.clear();
	getline(ifs, line);
	ifs.close();
	ofs.open("target_bf.txt", std::ios::app);
	ofs << line << std::endl;
	for (auto &m : map) {
		std::string out = bfdecoupling(&(m.second), &targetmol);
		if (out == "Yes") {
			ofs << m.first << std::endl;
		}
	}
	ofs.close();
}
