/*
 * rankscaffold.cpp
 *
 *  Created on: 2019年6月6日
 *      Author: yxchen
 */
/*
 * read substrate, coenzyme and certain protein; align substrate into existing
 * pocket by aligning common coenzyme molecule (e.g. P450) (or defined positions);
 * rank clashes for many possible rotations of the substrate so as to evaluate this
 * protein scaffold. (Clashes will not calculate sidechain farther than Cb)
 */
#include "basicfragment.h"
#include "mmmatches.h"
#include "atomtypessmarts.h"
#include "tmpltssas.h"
#include "cliquer/cliquer.h"
#include <iostream>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include "proteinrep/aaconformer.h"
#include "subsite.h"
using namespace myobcode;
using namespace subsitedesign;
using namespace OpenBabel;
using namespace NSPproteinrep;
#define DIST2CLASH 4.0

int main(int argc, char **argv)
{
	// read par
	std::string par(argv[1]);
	std::ifstream ifpar;
	ifpar.open(par.c_str());
	if (!ifpar.good())
	{
		std::cout << "no par" << std::endl;
		exit(1);
	}
	std::map<std::string, std::string> parmap;
	while(true)
	{
		std::string line;
		line.clear();
		getline(ifpar, line);
		if(!ifpar.good()) break;
		if (line.size()==0) continue;
		std::vector<std::string> words;
		std::stringstream input(line);
		std::string word;
		while(input>>word) words.push_back(word);
		parmap.insert(std::pair<std::string, std::string>(words[0], words[2]));
	}
	ifpar.close();

	const std::map<std::string, BasicFragment> &map =
			BasicFragment::getmap();
	OBConversion obconversion;
	OBMol tarcomol; // target coenzyme
	OBMol tarsubmol; // target substrate
	OBMol procomol; // protein coenzyme (and link to its protein scaffold by PDBid)
	obconversion.SetInFormat("sdf");
	std::string tarcosdf;
	std::string tarsubsdf;
	std::string procosdf;
	auto iter_par = parmap.find("target_coenzyme");
	if (iter_par != parmap.end())
		tarcosdf = iter_par->second;
	else
	{
		std::cout << "need target_coenzyme in par file" << std::endl;
		exit(1);
	}
	iter_par = parmap.find("target_substrate");
	if (iter_par != parmap.end())
		tarsubsdf = iter_par->second;
	else
	{
		std::cout << "need target_substrate in par file" << std::endl;
		exit(1);
	}
	iter_par = parmap.find("protein_coenzyme");
	if (iter_par != parmap.end())
		procosdf = iter_par->second;
	else
	{
		std::cout << "need proein_coenzyme in par file" << std::endl;
		exit(1);
	}
	// read the molecules
	obconversion.ReadFile(&tarcomol, tarcosdf);
	obconversion.ReadFile(&tarsubmol, tarsubsdf);
	obconversion.ReadFile(&procomol, procosdf);
	std::map<std::string, std::vector<std::shared_ptr<FMMatches>>> fragtarcomatches;
	std::map<std::string, std::vector<std::shared_ptr<FMMatches>>> fragprocomatches;
	for (auto &m : map) {
		fragtarcomatches[m.first] = findfmmatches(&(m.second), &tarcomol,
				false);
		fragprocomatches[m.first] = findfmmatches(&(m.second), &procomol,
					true);
	}
	std::vector<std::vector<std::string>> tarcoatomtypes_all=findatomtypes(&tarcomol);
	std::vector<std::string> tarcoatomtypes;
	for(auto &tat:tarcoatomtypes_all){
		if(tat.empty()) tarcoatomtypes.push_back("unspecified");
		else tarcoatomtypes.push_back(tat[0]);
	}
	std::vector<std::shared_ptr<SubstrAlignment>> ssas;
	ssas = findssas(fragprocomatches, fragtarcomatches);
	if (ssas.size() == 0)
	{
		std::cout << "no ssas" << std::endl;
		exit(1);
	}

	std::vector<std::string> mc_cb = {"N", "CA", "C", "O", "CB"};
	// read protein scaffold <- procomol.title and calculate compatibility
	std::vector<std::string> parsedname = TmpltSSAs::parsetmpltname(
			ssas[0]->mol1()->GetTitle()); // first line of procomol
	int modelno=std::atoi(parsedname[TmpltSSAs::MODELNO].c_str());
	char chainid=parsedname[TmpltSSAs::CHAINID][0];
	int residuenumber=std::atoi(parsedname[TmpltSSAs::RESIDUENO].c_str());
	auto reader=readpdbmodel(parsedname[TmpltSSAs::PDBID],1,modelno);

	if (reader == nullptr)
	{
		std::cout << "find an empty: " << parsedname[TmpltSSAs::PDBID] << std::endl;
		exit(1);
	}
/*
	// read scaffold from dir
	obconversion.SetInFormat("pdb");
	std::string scaffold;
	iter_par = parmap.find("scaffold");
	if (iter_par != parmap.end())
		scaffold = iter_par->second;
	else
	{
		std::cout << "need scaffold in par file" << std::endl;
		exit(1);
	}
	OBMol scaffoldmol; // scaffold
	obconversion.ReadFile(&scaffoldmol, scaffold);
*/
	obconversion.SetOutFormat("pdb");
	auto tarcossas=maketmpltssas(ssas,tarcoatomtypes);
	int ssaidx=0;
	std::map<int, std::vector<int>> clash_count; //[clash_count][move_id]
	for(TmpltSSAs::Alignment &algn:tarcossas->alignments){
		OBMol copy=tarsubmol;
		Move3D mv=algn.move3d;
		//mv.tvec[0]+=0.1;  ///<shift the aligned structure a little bit for visualization
		moveobmol(&copy, &mv);
		std::string outfile="Move_"+std::to_string(ssaidx++)+".pdb";
		obconversion.WriteFile(&copy, outfile);

		// fetch crd of copy (tar_sub)
		std::vector<NSPgeometry::XYZ> crds_ts;
		double *c = copy.GetConformer(0);
		for(int a = 0; a < copy.NumAtoms(); a++)
	    {
			NSPgeometry::XYZ cts(c[a*3], c[a*3+1], c[a*3+2]);
			crds_ts.push_back(cts);
	    }

		// calculate clash for this move
		int cc = 0;
		MapPdbKeyInt & mki = *(reader->mappdbkeyint());
		int ichain = mki.chainNumber(chainid);
		//reskey is pair of resiude sequence and insertion code from PDB record
		int iresidue = mki.posiNumber(std::pair<int,char>(residuenumber,' '), chainid);
		if (iresidue != -10)
		{
			AAConformersInModel cim;
			cim.getconformers(*reader);
			for (auto &chain : cim.conformers)
				for (auto &res : chain)
				{
					auto crd = res.getglobalcrd();
					for (auto &c : crd)
						if (find(mc_cb.begin(), mc_cb.end(), c.first) != mc_cb.end())
						// calculate clash
							for (auto &cts : crds_ts)
								if ((cts-c.second).squarednorm() < DIST2CLASH)
									cc++;
				}
		}
		else
			std::cout << "wrong resnum " << residuenumber << " in " << chainid << std::endl;
		if (clash_count.find(cc) != clash_count.end())
			clash_count[cc].push_back(ssaidx-1);
		else
		{
			std::vector<int> newkey = {ssaidx-1};
			clash_count.insert(std::pair<int, std::vector<int>>(cc, newkey));
		}
	}
	std::ofstream ofs("clash_rank.txt");
	ofs << "Move_id clash_num" << std::endl;
	for (auto &c : clash_count)
		for (auto cs : c.second)
			ofs << cs << "       " << c.first << std::endl;
}
