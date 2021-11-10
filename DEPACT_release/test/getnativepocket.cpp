/*
 * getnativepocket.cpp
 *
 *  Created on: 2019年5月5日
 *      Author: yxchen
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
#include "dataio/inputlines.h"
#include "dataio/datapaths.h"

/*
 * read .gnp , which contains all the infos about ligand, protein, wdmldpl and IsStrong
 * output the native pocket.
 */

using namespace myobcode;
using namespace subsitedesign;
using namespace OpenBabel;
using namespace NSPproteinrep;

static std::vector<std::string>  extractmodellines(const std::string &filename,int modelno){
	std::vector<std::string> res;
	NSPdataio::TextLines lines;
	lines.init(filename);
	bool inmodel = false;
	if (modelno == 0 ||modelno ==1)
		inmodel = true;
	int natoms = 0;
	for (auto &line : lines) {
		if (line.substr(0, 6) != "ATOM  " && line.substr(0, 6) != "HETATM"
				&& line.substr(0, 6) != "ENDMDL"
				&& line.substr(0, 5) != "MODEL") {
			res.push_back(line);
			continue;
		}
		if (line.substr(0, 5) == "MODEL") {
			inmodel = false;
			int num = std::stoi(line.substr(6));
			if (num == modelno)
				inmodel = true;
			continue;
		}
		if (line.substr(0, 6) == "ENDMDL") {
			inmodel = false;
			continue;
		}
		if (inmodel) {
			res.push_back( line);
			++natoms;
		}
	}
	return res;
}

int main(int argc, char **argv)
{

	// get contacts from native.pdb
	std::string ligand, protein;
	double wdmldpl = 1.0;
	bool IsWhole = true;
	std::ifstream ifpar(argv[1]); // gnp.gnp
	if(!ifpar.good()) {
		std::cout << ".gnp file failure" << std::endl;
		exit(1);
	}
	while(true) {
		std::string line;
		line.clear();
		getline(ifpar, line);
		if(!ifpar.good()) break;
		if (line.size()==0) continue;
		std::vector<std::string> words;
		std::stringstream input(line);
		std::string word;
		while(input>>word) words.push_back(word);
		if (words[0] == "Ligand" && words.size() == 3)
			ligand = words[2];
		if (words[0] == "Protein" && words.size() == 3)
			protein = words[2];
		if (words[0] == "Wdmldpl" && words.size() == 3)
			wdmldpl = std::stod(words[2]);
		if (words[0] == "IsStrongPocket" && words.size() == 3)
			if (words[2] == "True") IsWhole = false;
	}
	ifpar.close();

	const std::map<std::string, BasicFragment> &map =
			BasicFragment::getmap();
	OBConversion obconversion;
	obconversion.SetInFormat("sdf");
	std::ifstream ifs(ligand);
	obconversion.SetInStream(&ifs);
	OBMol nativemol;
	obconversion.Read(&nativemol);
	std::map<std::string, std::vector<std::shared_ptr<FMMatches>>> fragmatches;
	for (auto &m : map)
		fragmatches[m.first] = findfmmatches(&(m.second), &nativemol,
					true);
	std::vector<std::vector<std::string>> nativeatomtypes_all=findatomtypes(&nativemol);
	std::vector<std::string> nativeatomtypes;
	for(auto &nat:nativeatomtypes_all)
	{
		if(nat.empty()) //nativeatomtypes.push_back("unspecified");
		{
			std::cout << "Certain Atom(s) in the ligand cann't be decoded by atomtypesmarts.txt." << std::endl;
			exit(1);
		}
		else nativeatomtypes.push_back(nat[0]);
	}
	std::vector<std::shared_ptr<SubstrAlignment>> ssas;
	for(auto &m:map) {
		const std::vector<std::shared_ptr<FMMatches>> &fmmtarget=fragmatches.at(m.first);
		if(fmmtarget.empty()) continue;
		for (auto mtgt : fmmtarget) {
				auto nssa = std::shared_ptr < SubstrAlignment
				> (new SubstrAlignment(mtgt.get(), mtgt.get()));
				if (nssa->rmsd() > 1.0)
				continue;
				ssas.push_back(nssa);
		}
	}
// nativepocket = combine(ssas)
///*
	auto nssa = ssas[0];
		for (int j = 1; j < ssas.size(); j++)
			if (combinable(*nssa, *(ssas[j])))
				nssa = std::shared_ptr<SubstrAlignment> (new SubstrAlignment(*nssa, *(ssas[j])));
	ssas.clear();
	ssas.push_back(nssa);
//*/

	auto tmpltssas=maketmpltssas(ssas,nativeatomtypes);
	std::vector<std::string> parsedname = TmpltSSAs::parsetmpltname(
			tmpltssas->tmpltname);
	int modelno=std::atoi(parsedname[TmpltSSAs::MODELNO].c_str());
	char chainid=parsedname[TmpltSSAs::CHAINID][0];
	int residuenumber=std::atoi(parsedname[TmpltSSAs::RESIDUENO].c_str());

	std::shared_ptr<PdbReader> reader=std::shared_ptr<PdbReader> (new PdbReader());
	if (!protein.empty())
	{
		std::vector<std::string> lines=extractmodellines(protein,modelno);
		reader->readpdb(lines);
	}
	else
		reader=readpdbmodel(parsedname[TmpltSSAs::PDBID],1,modelno);

	MapPdbKeyInt & mki=*(reader->mappdbkeyint());
	int ichain=mki.chainNumber(chainid);
	//reskey is pair of resiude sequence and insertion code from PDB record
	int iresidue=mki.posiNumber(std::pair<int,char>(residuenumber,' '),chainid);
	if (iresidue == -10)
	{
		std::cout << "this pdb cannot be read." << std::endl;
		exit(1);
	}
	AAConformersInModel cim;
	cim.getconformers(*reader);
	auto contacts=findcontacts(&cim,ichain,iresidue);
	// make native_pocket
	std::ostringstream pdbsstr;
	obconversion.SetOutFormat("pdb");
	obconversion.Write(&nativemol, &pdbsstr);
	std::istringstream isstrpdb(pdbsstr.str());
	TargetStruct tgt(nativemol.GetTitle(), nativeatomtypes, isstrpdb);
	std::vector<std::shared_ptr<Subsite>> nativessites;

	int q = 0;
	for (auto &aln : tmpltssas->alignments)
	{
		std::shared_ptr<Subsite> ss(new Subsite(tgt,*contacts,aln, IsWhole));
		if (ss->keyresidues().empty()) continue;
		nativessites.push_back(ss);
	}

	double min_s = 0.0;
	int idx = 0;
	std::string outns = "nativescore.txt";
	std::ofstream ofnatscore(outns.c_str(), std::ios::app);
	for (int i = 0; i < nativessites.size(); i++)
	{
		auto &nsub = nativessites[i];
		std::string nfile = "Native_"+std::to_string(i)+".pdb";
		nsub->writepdb(nfile);
		std::vector<double> scores = nsub->detailtotalscore(outns, wdmldpl);
		if (min_s > scores[0])
		{
			min_s = scores[0];
			idx = i;
		}
		ofnatscore << "Native " << std::to_string(i) << " total score: " << scores[0] << std::endl;
	}
	ofnatscore << "the best pocket for native: " << idx << std::endl;
}


