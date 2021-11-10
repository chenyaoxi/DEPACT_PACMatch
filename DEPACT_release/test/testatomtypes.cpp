/*
 * testatomtypes.cpp
 *
 *  Created on: 2018年8月16日
 *      Author: hyliu
 */

#include "atomtypessmarts.h"
using namespace myobcode;
using namespace OpenBabel;
/**
 * @brief test whether the atom type definitions can cover ligands read from input
 *
 * OpenBabel tools are used to input/output ligand and to find SMARTS matches
 * The first argument is the filename of sdf contains the input ligands
 * The SMARTS patterns defining atomtypes are read from "atomtypesmarts.txt"
 */
int main(int argc, char **argv) {
	OBConversion obconversion;  ///<OpenBabel object to read andconvert molecule
	OBConversion conv1;
	std::ofstream ofs("unknownatoms.txt");  ///<output uncovered ligand atoms
	obconversion.SetInFormat("sdf");
	obconversion.SetOutFormat("pdb");
	conv1.SetOutFormat("sdf");
	std::string molsdf(argv[1]);
	std::ifstream ifs(molsdf);
	obconversion.SetInStream(&ifs);
	while(ifs){
		OBMol mol;
		obconversion.Read(&mol);
		std::vector<std::vector<std::string>> atypes=findatomtypes(&mol);
		int idx=0;
		std::cout <<mol.GetTitle()<<std::endl;
		bool ok=true;
		for(auto &s:atypes){
			std::cout <<"atom "<<++idx;
			for(auto &t:s)std::cout <<" "<<t;
			std::cout<<std::endl;
			if(s.empty() || s.size()> 1) {
				ok=false;
			}
		}
		if(!ok){
			int nhvyatoms=mol.NumHvyAtoms();
			if(nhvyatoms<4) continue;
			std::string file=std::string(mol.GetTitle()).substr(0,8)+".pdb";
			std::string filesdf=std::string(mol.GetTitle()).substr(0,8)+".sdf";
			obconversion.WriteFile(&mol,file);
			conv1.WriteFile(&mol,filesdf);
			idx=0;
			for(auto &s:atypes){
					ofs<<mol.GetTitle() <<" atom "<<++idx;
					if(s.empty()) ofs <<"unspecified";
					for(auto &t:s){
						ofs<<" "<<t;
					}
					ofs<<std::endl;
			}
		}
	}
}
