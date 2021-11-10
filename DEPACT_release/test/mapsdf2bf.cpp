/*
 * mapsdf2bf.cpp
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
//all-sdf.sdf
	const std::map<std::string, BasicFragment> &map =
			BasicFragment::getmap();

	std::string sdf(argv[1]);
	std::ifstream ifs;
	ifs.open(sdf.c_str());
	if (!ifs.good()) {
		std::cout << "sdf failure" << std::endl;
		exit(1);
	}
	int start = 1;
	std::string single_sdf = "single.sdf";
	std::ofstream ofs;
	std::string name;
	while(true) {
		std::string line;
		line.clear();
		getline(ifs, line);
		if(!ifs.good()) break;
		if (line.size()==0) continue;
		if (start == 1)
			ofs.open(single_sdf.c_str());
		else
			ofs.open(single_sdf.c_str(),std::ios::app);
		ofs << line << std::endl;
		ofs.close();
		if (line != "$$$$") {
			if (start == 1) {
				name = line;
				start = 0;
			}
		}
		else {
			start = 1;
			OBConversion obconversion;
			OBMol mol;
			obconversion.SetInFormat("sdf");
			obconversion.ReadFile(&mol, single_sdf);
			std::map<std::string, std::string> map_sdf2bf;
			for (auto &m : map) {
				std::string out = bfdecoupling(&(m.second), &mol);
				if (out == "Yes") {
					std::string bf = "bf_" + m.first + ".txt";
					std::ofstream ofbf;
					ofbf.open(bf.c_str(), std::ios::app);
					ofbf << name << std::endl;
					ofbf.close();
				}
			}
		}
	}
	ifs.close();
//should write bf-sdf.txt

}
