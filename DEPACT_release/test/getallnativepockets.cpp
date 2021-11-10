/*
 * getallnativepockets.cpp
 *
 *  Created on: 2020年7月28日
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
 * input: list.txt
 * extract all native pockets in the list.txt
 * output: Native_i.pdb nativescore.txt(also store the list_line to know PDBID_ligand name)
 */

using namespace std;
using namespace myobcode;
using namespace subsitedesign;
using namespace OpenBabel;
using namespace NSPproteinrep;
using namespace NSPdataio;

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
	// prepare
	const std::map<std::string, BasicFragment> &map =
			BasicFragment::getmap();
	OBConversion obconversion;
	obconversion.SetInFormat("sdf");

	// read every line from line.txt
	std::ifstream ifl(argv[1]); // line.txt
	vector<string> lines;
	if (!ifl.good())
	{
		std::cout << "list.txt is wrong" << std::endl;
		exit(1);
	}
	while (true)
	{
		string line;
		line.clear();
		getline(ifl, line);
		if (!ifl.good()) break;
		if (line.size() == 0) continue;
		lines.push_back(line);
	}

	bool all_include;
	string input2 = argv[2];
	if (input2 == "True")
		all_include = true;
	else if (input2 == "False")
		all_include = false;
	else
	{
		std::cout << "second parameter should be True or False" << std::endl;
		exit(1);
	}


	string as = getenvpath("DEPACT_DATAPATH")+"all-sdf.sdf";
	std::ifstream ifa(as.c_str()); // all-sdf.sdf
	if (!ifa.good())
	{
		std::cout << "all-sdf.sdf is wrong" << std::endl;
		exit(1);
	}
	bool find = false;
	bool start = true;
	std::ofstream oft;
	int ni = 0; // Native_i.pdb
	string recent_name;
	std::string outns = "native_detail.txt";
	std::ofstream ofnatscore(outns.c_str(), std::ios::app);
	std::string outsc = "native_score.txt";
	std::ofstream ofnsc(outsc.c_str(), std::ios::app);
	ofnsc << "fullname   pocket_id   score" << std::endl;
	std::pair<double, int> min_poc; // scorr, id
	min_poc.first = 0.0;
	min_poc.second = -1;
	double wdmldpl = 1.0;
	if (argc < 4)
		std::cout << "wdmldpl = 1.0, otherwise the third input parameter as weight is required" << std::endl;
	else wdmldpl = std::stod(argv[3]);
	while (true)
	{
		std::string line;
		line.clear();
		getline(ifa, line);
		if (!ifa.good()) break;
		if (line.size() == 0) continue;
		for (int i = 0; i < lines.size(); i++)
			if (lines[i] == line)
			{
				find = true;
				recent_name = line;
			}
		if (find)
		{
			if (start)
				oft.open("tarnat.sdf");
			else
				oft.open("tarnat.sdf", std::ios::app);
			oft << line << std::endl;
			oft.close();
			if (line != "$$$$")
			{
				if (start)
					start = false;
			}
			else
			{
				start = true;
				find = false;
				OBConversion obconversion;
				OBMol nativemol;
				obconversion.SetInFormat("sdf");
				obconversion.ReadFile(&nativemol, "tarnat.sdf");
				std::map<std::string, std::vector<std::shared_ptr<FMMatches>>> fragmatches;
				for (auto &m : map)
					fragmatches[m.first] = findfmmatches(&(m.second), &nativemol,
								true);
				std::vector<std::vector<std::string>> nativeatomtypes_all=findatomtypes(&nativemol);
				std::vector<std::string> nativeatomtypes;
				for(auto &nat:nativeatomtypes_all)
				{
					if(nat.empty()) nativeatomtypes.push_back("unspecified");
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
				// ssas can be empty.
				if (ssas.size() == 0) continue;
				auto nssa = ssas[0];
					for (int j = 1; j < ssas.size(); j++)
						if (combinable(*nssa, *(ssas[j])))
							nssa = std::shared_ptr<SubstrAlignment> (new SubstrAlignment(*nssa, *(ssas[j])));
				ssas.clear();
				ssas.push_back(nssa);
				auto tmpltssas=maketmpltssas(ssas,nativeatomtypes);
				std::vector<std::string> parsedname = TmpltSSAs::parsetmpltname(
						tmpltssas->tmpltname);
				int modelno=std::atoi(parsedname[TmpltSSAs::MODELNO].c_str());
				char chainid=parsedname[TmpltSSAs::CHAINID][0];
				int residuenumber=std::atoi(parsedname[TmpltSSAs::RESIDUENO].c_str());

				std::shared_ptr<PdbReader> reader=std::shared_ptr<PdbReader> (new PdbReader());
				reader=readpdbmodel(parsedname[TmpltSSAs::PDBID],1,modelno);
				if (reader == nullptr)
					std::cout << "find an empty: " << parsedname[TmpltSSAs::PDBID] << std::endl;
				else
				{
					MapPdbKeyInt & mki=*(reader->mappdbkeyint());
					int ichain=mki.chainNumber(chainid);
					//reskey is pair of resiude sequence and insertion code from PDB record
					int iresidue=mki.posiNumber(std::pair<int,char>(residuenumber,' '),chainid);
					if (iresidue == -10) continue;
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
	//				vector<string> codes;
					int q = 0;
	//				int tgt_natm = tgt.conformer.atomlist.size();
					for (auto &aln : tmpltssas->alignments)
					{
						std::shared_ptr<Subsite> ss(new Subsite(tgt,*contacts,aln, all_include));
						if (ss->keyresidues().empty()) continue;
						nativessites.push_back(ss);
						// getting code
	//					int codenum = 0;
	//					string code(tgt_natm, '0');
	//					for (auto ap : aln.alignedpairs)
	//					{
	//						code[ap.first] = '1';
	//						codenum--;
	//					}
	//					codes.push_back(code);
					}
					if (nativessites.size() > 0)
					{
						double min_s = 1000000000.0;
						int idx = 0;
						for (int i = 0; i < nativessites.size(); i++)
						{
							auto &nsub = nativessites[i];
							double score = nsub->totalscore(wdmldpl);
							// change...
							std::ofstream ofs(outns.c_str(), std::ios::app);
							ofs << "Native_"+std::to_string(ni)+".pdb" << endl;
							std::vector<double> scores = nsub->detailtotalscore(outns, wdmldpl);
							score = scores[0];
							if (min_s > score)
							{
								min_s = score;
								idx = i;
							}
						}
						std::string nfile = "Native_"+std::to_string(ni)+".pdb";
						auto nsub = nativessites[idx];
						nsub->writepdb(nfile);
//						ofnatscore << "Native_" << std::to_string(ni) << " total score " << min_s << std::endl;
		//				int count_c = 0;
		//				for (auto c : codes[idx])
		//					if (c == '1')
		//						count_c++;
						ofnsc << recent_name << " Native_" << std::to_string(ni) << ".pdb " << min_s << std::endl;
						if (min_s < min_poc.first)
						{
							min_poc.first = min_s;
							min_poc.second = ni;
						}
						ni++;
					} // native have contacts.
				} // read successfully
			} // start getting native pocket
		} // find a native protein
	} // read all-sdf.sdf
	ifa.close();
	ofnsc << "min_score pocket is Native_" << std::to_string(min_poc.second)
	    << ".pdb with score " << min_poc.first << std::endl;
	std::cout << "finish all process" << std::endl;


}




