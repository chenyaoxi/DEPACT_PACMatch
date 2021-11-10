/*
 * testsmartsfragment.cpp
 *
 *  Created on: 2018年8月6日
 *      Author: hyliu
 */
#include "basicfragment.h"
#include "mmmatches.h"
#include "atomtypessmarts.h"
#include "tmpltssas.h"
#include "cliquer/cliquer.h"
#include <iostream>
#include "buildpocket.h"

//#include "../openbabel/include/openbabel/mol.h"
//#include "../openbabel/include/openbabel/obconversion.h"
using namespace myobcode;
using namespace subsitedesign;
using namespace OpenBabel;
static double rmsdmax(int natoms){
	if(natoms<=4) return 0.1;
	else if(natoms <=8) return 0.25;
//	else if(natoms <=16) return 0.4;
	else return 0.5;
}
void printalignment(const SubstrAlignment &ssa) {
	for (int m = 0; m < ssa.size(); ++m) {
		auto pair = ssa.alignedpair(m);
//		if (ssa.atomspecial(m))
//			std::cout << " S";
//		else
			std::cout << "  ";
		std::cout << pair.first << ":" << pair.second;
	}
	std::cout << std::endl;
}
/**
 * test substructure alignment codes between two ligand molecules
 *
 * argv[1] is the sdf file containing the target ligand
 * argv[2] is the sdf file containing the template
 * basic fragments and atom types defined using SMARTS are read from
 * "basicfragments.txt" and "atomtypesmarts.txt", respectively.
 *
 * The substructure alignments are written to standard output in XML format,
 * for later use @see testxmlio.cpp
 * For verification purpose, the template are superimposed with target
 * according to each substructure alignments and the results written out in PDB format
 */

int main(int argc, char **argv) {
	const std::map<std::string, BasicFragment> &map =
			BasicFragment::getmap();
	/*	for(auto &m:map){
	 std::cout <<m.first <<" "<<m.second.smarts;
	 for(auto i:m.second.scatoms) std::cout <<' '<<i;
	 std::cout<<std::endl;
	 }*/
	OBConversion obconversion;
	OBMol targetmol;
	OBMol templatemol;
	obconversion.SetInFormat("sdf");
	std::string targetsdf(argv[1]);
	std::string templatesdf(argv[2]);
//	std::string outfile(argv[3]);
	///read the molecules using OBConversion
	obconversion.ReadFile(&targetmol, targetsdf);
	obconversion.ReadFile(&templatemol, templatesdf);
	std::map<std::string, std::vector<std::shared_ptr<FMMatches>>> fragtargetmatches;
	std::map<std::string, std::vector<std::shared_ptr<FMMatches>>> fragtmpltmatches;
	///find all matches to all basic fragments
	for (auto &m : map) {
		fragtargetmatches[m.first] = findfmmatches(&(m.second), &targetmol,
				false);
		fragtmpltmatches[m.first] = findfmmatches(&(m.second), &templatemol,
					true);
	}
	///determine atom types of target ligand
	std::vector<std::vector<std::string>> targetatomtypes_all=findatomtypes(&targetmol);
	std::vector<std::string> targetatomtypes;
	for(auto &tat:targetatomtypes_all){
		if(tat.empty()) targetatomtypes.push_back("unspecified");
		else targetatomtypes.push_back(tat[0]);
	}
	///prepare output XML DOM
//	std::string initdoc("<root/>");
//	xmlpp::DomParser parser;
//	parser.parse_memory(initdoc);
//	xmlpp::Node * root=parser.get_document()->get_root_node();
//	addnodetargetstruct(root,&targetmol,targetatomtypes);
	std::ostringstream pdbsstr;
	obconversion.SetOutFormat("pdb");
	obconversion.Write(&targetmol,&pdbsstr);
	std::istringstream isstrpdb(pdbsstr.str());
	TargetStruct tgs(targetmol.GetTitle(),targetatomtypes,isstrpdb);
	int atmsize = tgs.conformer.atomlist.size();
//	std::ofstream ofs(outfile);
//	{
//		boost::archive::text_oarchive oa(ofs);
//		oa << tgs;
//	}
	///find substructure alignments
	std::vector<std::shared_ptr<SubstrAlignment>> ssas;
//	ssas = findssas(fragtargetmatches, fragtmpltmatches);
	std::vector<std::string> smartsid;
	for(auto &m:map) {
		const std::vector<std::shared_ptr<FMMatches>> &fmmtarget=fragtargetmatches.at(m.first);
		if(fmmtarget.empty()) continue;
		const std::vector<std::shared_ptr<FMMatches>> &fmmtemplate =
		fragtmpltmatches.at(m.first);
		for (auto mtgt : fmmtarget) {
			for (auto mtmplt : fmmtemplate) {
				auto nssa = std::shared_ptr < SubstrAlignment
				> (new SubstrAlignment(mtgt.get(), mtmplt.get()));
				if (nssa->rmsd() > rmsdmax(nssa->size()))
				continue;
				/*			int idx = 0;
				 bool redundant = false;
				 for (auto &ossa : ssas) {
				 if (equivalent(*nssa, *ossa)) {
				 ossa = std::shared_ptr < SubstrAlignment
				 > (new SubstrAlignment(*nssa, *ossa));
				 smartsid[idx] = smartsid[idx] + "_M_" + m.first;
				 redundant = true;
				 break;
				 }
				 idx++;
				 }
				 if (!redundant) {*/
				ssas.push_back(nssa);
				smartsid.push_back(m.first);
//					idx++;}
			}
		}
	}

	// output bf_ssa_relationship.
	std::ofstream ofb("bf_ssa.txt");
	for (int i = 0; i < smartsid.size(); i++)
		ofb << i << " " << smartsid[i] << std::endl;

	///determine atomtypes of template
	std::vector<std::vector<std::string>> tmpltatomtypes_all=findatomtypes(&templatemol);
	std::vector<std::string> tmpltatomtypes;
	for(auto &tat:tmpltatomtypes_all){
		if(tat.empty()) tmpltatomtypes.push_back("unspecified");
		else tmpltatomtypes.push_back(tat[0]);
	}

	///wrap up the substruature alignments and other data
	auto tmpltssas=maketmpltssas(ssas,tmpltatomtypes);
	int ntmplts=1;
//	ofs << ntmplts <<std::endl; //number of objects
//	{
//		boost::archive::text_oarchive oa(ofs);
//		oa <<*tmpltssas;
//	}
	///write XML node
//	addnodetmpltssas(root,*tmpltssas);
//	parser.get_document()->write_to_stream_formatted(std::cout);

	// check if all atoms are covered.
	std::vector<bool> atmscover(atmsize, false);
	///write the transformed templates in PDB format
	obconversion.SetOutFormat("pdb");
	int ssaidx=0;
	std::ofstream oftxt("alignment.txt");
	for(TmpltSSAs::Alignment &algn:tmpltssas->alignments){
		oftxt << "SSA_" << ssaidx << ":" << std::endl;
		for (auto a : algn.alignedpairs)
		{
			oftxt << " " << a.first << " " << a.second << std::endl;
			atmscover[a.first] = true;
		}
		OBMol copy=templatemol;
		Move3D mv=algn.move3d;
		//mv.tvec[0]+=0.1;  ///<shift the aligned structure a little bit for visualization
		moveobmol(&copy, &mv);
		std::string outfile="SSA_"+std::to_string(ssaidx++)+".pdb";
		obconversion.WriteFile(&copy, outfile);
	}
	bool allcover = true;
	for (int a = 0; a < atmscover.size(); a++)
		if (!atmscover[a])
		{
			oftxt << a << " is not covered" << std::endl;
			allcover = false;
		}
	if (allcover)
		oftxt << "all atoms are covered" << std::endl;
	oftxt.close();
}
