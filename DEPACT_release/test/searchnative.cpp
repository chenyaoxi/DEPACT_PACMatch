/*
 * searchnative.cpp
 *
 *  Created on: 2019年4月15日
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

void makessites(TmpltSSAs &ssa, TargetStruct &tgt,
		std::vector<std::shared_ptr<Subsite>> &ssites)
{
	std::vector<std::string> parsedname = TmpltSSAs::parsetmpltname(
			ssa.tmpltname);
	int modelno=std::atoi(parsedname[TmpltSSAs::MODELNO].c_str());
	char chainid=parsedname[TmpltSSAs::CHAINID][0];
	int residuenumber=std::atoi(parsedname[TmpltSSAs::RESIDUENO].c_str());
	auto reader=readpdbmodel(parsedname[TmpltSSAs::PDBID],1,modelno);
	if (reader == nullptr)
		std::cout << "find an empty: " << parsedname[TmpltSSAs::PDBID] << std::endl;
	else
	{
		MapPdbKeyInt & mki = *(reader->mappdbkeyint());
		int ichain = mki.chainNumber(chainid);
		//reskey is pair of resiude sequence and insertion code from PDB record
		int iresidue = mki.posiNumber(std::pair<int,char>(residuenumber,' '), chainid);
		if (iresidue != -10)
		{
			AAConformersInModel cim;
			cim.getconformers(*reader);
			auto contacts = findcontacts(&cim,ichain,iresidue);
			for (auto &aln : ssa.alignments)
			{
				std::shared_ptr<Subsite> ss(new Subsite(tgt,*contacts,aln));
				if (ss->keyresidues().empty()) continue;
				ssites.push_back(ss);
			}
		}
		else
			std::cout << "wrong resnum " << residuenumber << " in " << chainid << std::endl;
	}
}

int main(int argc, char **argv)
{
//ligand_name (ATP)
//read target and tmplts
	const std::map<std::string, BasicFragment> &map =
				BasicFragment::getmap();
	OBConversion obconversion;
	OBMol targetmol;
	OBMol templatemol;
	obconversion.SetInFormat("sdf");
	std::string ligand(argv[1]);
	obconversion.ReadFile(&targetmol, ligand+".sdf"); // ATP.sdf
	std::map<std::string, std::vector<std::shared_ptr<FMMatches>>> fragtargetmatches;
	for (auto &m : map)
		fragtargetmatches[m.first] = findfmmatches(&(m.second), &targetmol,
				false);
	///determine atom types of target ligand
	std::vector<std::vector<std::string>> targetatomtypes_all=findatomtypes(&targetmol);
	std::vector<std::string> targetatomtypes;
	for(auto &tat:targetatomtypes_all)
	{
		if(tat.empty()) targetatomtypes.push_back("unspecified");
		else targetatomtypes.push_back(tat[0]);
	}
	std::ostringstream pdbsstr;
	obconversion.SetOutFormat("pdb");
	obconversion.Write(&targetmol,&pdbsstr);
	std::istringstream isstrpdb(pdbsstr.str());
	TargetStruct tgt(targetmol.GetTitle(),targetatomtypes,isstrpdb);

	std::string t_file("tmplt.sdf");
	std::ifstream ifall;
	ifall.open("all-sdf.sdf");
	int find = 0;
	std::ofstream of_np("native_pdb.txt", std::ios::app);
	while (true)
	{
		std::string line;
		line.clear();
		getline(ifall, line);
		if(!ifall.good()) break;
		if (line.size() == 0) continue;
		if (line.find(ligand) == 5)
		{
			find = 1;
			of_np << line << std::endl;
		}
		if (find == 1)
		{
			std::ofstream of_t;
			of_t.open(t_file.c_str(), std::ios::app);
			of_t << line << std::endl;
			of_t.close();
		}
		if (find == 1 && line == "$$$$")
		{
			find = 0;
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
			for(auto &tat:tmpltatomtypes_all)
			{
				if(tat.empty()) tmpltatomtypes.push_back("unspecified");
				else tmpltatomtypes.push_back(tat[0]);
			}
			///wrap up the substruature alignments and other data
			if (ssas.size() == 0)
			{
				of_np << "ssas.size() = 0" << std::endl;
				remove(t_file.c_str());
				continue;
			}
			auto tmpltssa = maketmpltssas(ssas,tmpltatomtypes);
			//make ssites
			int ss_idx = 0;
			std::vector<std::shared_ptr<Subsite>> ssites;
			makessites(*tmpltssa, tgt, ssites);
			double sum_pl = 0.0;
			double sum_ml = 0.0;
			for (auto &ss : ssites)
			{
				for (auto & pl : ss->keyresiduescores().plscores)
					sum_pl += ss->keyresiduescores().totalplscore(pl.first);
				for (auto & ml : ss->keyresiduescores().mlscores)
					sum_ml += ss->keyresiduescores().totalplscore(ml.first);
			}
			of_np << "sum_pl: " << sum_pl << std::endl;
			of_np << "sum_ml: " << sum_ml << std::endl;
			remove(t_file.c_str());
		}
	}

	return 0;
}


