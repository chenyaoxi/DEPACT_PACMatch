/*
 * overall_subpool.cpp
 *
 *  Created on: 2019年4月1日
 *      Author: yxchen
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>
#include "basicfragment.h"
#include "mmmatches.h"
#include "atomtypessmarts.h"
#include "tmpltssas.h"
#include "cliquer/cliquer.h"
#include "subsite.h"
#include "proteinrep/aaconformer.h"

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
using namespace myobcode;
using namespace subsitedesign;
using namespace OpenBabel;
using namespace NSPproteinrep;

int main(int argc, char **argv) {
//target_bf.txt bf_sum target.sdf -> ssas.dat
//read target_bf.txt
	std::string target(argv[1]);
	std::ifstream ifs;
	ifs.open(target.c_str());
	if (!ifs.good()) {
		std::cout << "sdf failure" << std::endl;
		exit(1);
	}
	int start = 1;
	std::string target_name;
	std::vector<std::string> bflists;
	while(true) {
		std::string line;
		line.clear();
		getline(ifs, line);
		if(!ifs.good()) break;
		if (line.size()==0) continue;
		if (start == 1) {
			target_name = line;
			start = 0;
		}
		else
			bflists.push_back(line);
	}
	ifs.close();

//make map for template->bf_count
	std::map<std::string, int> bf_map;
	for (auto bf : bflists) {
		std::string bf_fn = "bf_"+bf+".txt";
		std::ifstream if_bf;
		if_bf.open(bf_fn.c_str());
		while(true) {
			std::string line;
			line.clear();
			getline(if_bf, line);
			if(!if_bf.good()) break;
			if (line.size()==0) continue;
			std::map<std::string, int>::iterator iter = bf_map.find(line);
			if (iter != bf_map.end())
				bf_map[line] = bf_map[line]+1;
			else
				bf_map.insert(std::pair<std::string, int>(line, 1));
		}
		if_bf.close();
	}
// if bf_sum > 1
	int bf_sum = std::stoi(argv[2]);
	std::vector<std::string> tmplts;
	for (auto &m : bf_map)
		if (m.second >= bf_sum)
			tmplts.push_back(m.first);
	std::cout << tmplts.size() << " templates will be analyzed" << std::endl;
// if bf_sum = 1: std::cout << bf_map.size() << " templates will be analyzed" << std::endl;

//read target and tmplts
	const std::map<std::string, BasicFragment> &map =
				BasicFragment::getmap();
	OBConversion obconversion;
	OBMol targetmol;
	OBMol templatemol;
	obconversion.SetInFormat("sdf");
	std::string targetsdf(argv[3]);
	obconversion.ReadFile(&targetmol, targetsdf);
	std::map<std::string, std::vector<std::shared_ptr<FMMatches>>> fragtargetmatches;
	for (auto &m : map)
		fragtargetmatches[m.first] = findfmmatches(&(m.second), &targetmol,
				false);
	///determine atom types of target ligand
	std::vector<std::vector<std::string>> targetatomtypes_all=findatomtypes(&targetmol);
	std::vector<std::string> targetatomtypes;
	for(auto &tat:targetatomtypes_all){
		if(tat.empty()) targetatomtypes.push_back("unspecified");
		else targetatomtypes.push_back(tat[0]);
	}
	std::ostringstream pdbsstr;
	obconversion.SetOutFormat("pdb");
	obconversion.Write(&targetmol,&pdbsstr);
	std::istringstream isstrpdb(pdbsstr.str());
	TargetStruct tgt(targetmol.GetTitle(),targetatomtypes,isstrpdb);
	std::ofstream ofs("ssas.dat");
	{
		boost::archive::text_oarchive oa(ofs);
		oa << tgt;
	}

	std::string t_file = "template.sdf";
//	std::vector<TmpltSSAs> tmpltssas;

	ofs << tmplts.size() << std::endl; // tmplts.size()
//establish tmpltssas
	for (int i = 0; i < tmplts.size(); i++) {
	//test ->ssites.txt
//	for (int i = 0; i < 10; i++) {
		//output single tmplt
		int find = 0;
		std::string temp = tmplts[i];
		remove("template.sdf");
		std::ifstream ifall;
		ifall.open("all-sdf.sdf");
		while(true) {
			std::string line;
			line.clear();
			getline(ifall, line);
			if(!ifall.good()) break;
			if (line.size() == 0) continue;
			if (line == temp) find = 1;
			if (find == 1) {
				std::ofstream of_t;
				of_t.open(t_file.c_str(), std::ios::app);
				of_t << line << std::endl;
			}
			if (find == 1 && line == "$$$$") break;
		}
		ifall.close();
		//for tmplt part
		obconversion.ReadFile(&templatemol, t_file);
		std::map<std::string, std::vector<std::shared_ptr<FMMatches>>> fragtmpltmatches;
		for (auto &m : map)
			fragtmpltmatches[m.first] = findfmmatches(&(m.second), &templatemol,
						true);
		///find substructure alignments
		std::vector<std::shared_ptr<SubstrAlignment>> ssas;
		ssas = findssas(fragtargetmatches, fragtmpltmatches);
		///determine atomtypes of template
		std::vector<std::vector<std::string>> tmpltatomtypes_all=findatomtypes(&templatemol);
		std::vector<std::string> tmpltatomtypes;
		for(auto &tat:tmpltatomtypes_all){
			if(tat.empty()) tmpltatomtypes.push_back("unspecified");
			else tmpltatomtypes.push_back(tat[0]);
		}
		///wrap up the substruature alignments and other data
		if (ssas.size() == 0) continue;
//		std::cout << ssas[0]->mol2()->GetTitle() << std::endl;
		auto tmpltssa = maketmpltssas(ssas,tmpltatomtypes);
//		tmpltssas.push_back(*tmpltssa);
		{
			boost::archive::text_oarchive oa(ofs);
			oa <<*tmpltssa;
		}
	}

}
