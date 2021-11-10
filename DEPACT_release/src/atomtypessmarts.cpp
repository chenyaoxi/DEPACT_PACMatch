/*
 * atomtypessmarts.cpp
 *
 *  Created on: 2018年8月16日
 *      Author: hyliu
 */
#include "atomtypessmarts.h"
#include "geometry/calculators.h"
using namespace subsitedesign;
using namespace OpenBabel;
std::vector<std::vector<std::string>> myobcode::findatomtypes(OpenBabel::OBMol *mol){
	    static std::vector<std::string> othertypes{"E_other","X_eother"};
		int natoms=mol->NumAtoms();
		std::vector<std::vector<std::string>> types(natoms,std::vector<std::string>());
		auto & map=AtomType::getmap();
		for(auto &at:map){
			OpenBabel::OBSmartsPattern smarts;
			if(at.first=="E_other" || at.first=="X_eother") continue;
			smarts.Init(at.second.smarts);
			std::vector<std::vector<int> > maplist;
			if (smarts.Match(*mol)) {
				maplist = smarts.GetUMapList();
				for (auto itr = maplist.begin(); itr != maplist.end(); itr++) {
					types[(*itr)[0]-1].push_back(at.second.name);
				}
			}
		}
		for(auto &t:othertypes){
			OpenBabel::OBSmartsPattern smarts;
			smarts.Init(map.at(t).smarts);
			std::vector<std::vector<int> > maplist;
			if (smarts.Match(*mol)) {
				maplist = smarts.GetUMapList();
				for (auto itr = maplist.begin(); itr != maplist.end(); itr++) {
					if(types[(*itr)[0]-1].empty())
						types[(*itr)[0]-1].push_back(t);
				}
			}
		}
		return types;
}
std::vector<std::vector<int>> myobcode::connectivity(OBMol *mol){
	std::vector<NSPgeometry::XYZ> crds;
	OBAtom *atom;
	int natoms=mol->NumAtoms();
	for(int i=0;i<natoms;++i){
		atom=mol->GetAtom(i+1);
		crds.push_back(NSPgeometry::XYZ(atom->x(),atom->y(),atom->z()));
	}
	std::vector<std::vector<int>> res(natoms,std::vector<int>());
	for(int i=0;i<natoms-1;++i){
		for(int j=i+1;j<natoms;++j){
			double rij2=(crds[i]-crds[j]).squarednorm();
			if(rij2<3.61){
				res[i].push_back(j);
				res[j].push_back(i);
			}
		}
	}
	return res;
}


