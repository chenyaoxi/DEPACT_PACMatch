/*
 * statis_metal_pdbs.cpp
 *
 *  Created on: 2019年12月19日
 *      Author: yxchen
 */

#include "buildpocket.h"
#include "statis_metal.h"
#include "proteinrep/pdbreader.h"
#include "proteinrep/aaconformer.h"
using namespace NSPproteinrep;
using namespace buildpocket;
using namespace statismetal;

// getting contact conditions (bond & angle) of metal from many_pdbs in one file.
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

// EveryPdb->output
	// from titles.txt, read every sdf's complex.pdb
	std::string pdbtitle = getp(parmap, "pdbs_titles");
	int num_th = std::stoi(getp(parmap, "num_th"));
	std::vector<std::string> pdbs;
	std::ifstream ifs;
	ifs.open(pdbtitle.c_str());
	if (!ifs.good())
	{
		std::cout << "no pdbs_title" << std::endl;
		exit(1);
	}
	while (true)
	{
		std::string line;
		getline(ifs, line);
		if (!ifs.good()) break;
		if (line.size() == 0) continue;
		pdbs.push_back(line);
	}
	ifs.close();
	std::string spec_metal = getp(parmap, "specific_metal");
	if (spec_metal == "" || spec_metal == "ALL")
		std::cout << "Statistics for all coordinating contacts." << std::endl;
	else
		std::cout << "Statistics for " << spec_metal << std::endl;

	int id = 0;
	for (auto pdb : pdbs)
	{
		id++;
		if (id > num_th)
			break;
		auto cim = getconfsfromtitlefile_spdb(pdb);
		// statistic: bond_length, angle and dihedral angle.
		if (spec_metal == "" || spec_metal == "ALL")
		{
			// extract metal and other atoms
			std::vector<MetalConf> mcs = aaconf2metalconf_single(cim);
			metalconf2statis_print(mcs);
		}
		else
		{
			// extract metal and other atoms
			std::vector<MetalConf> mcs = aaconf2metalconf_spec_single(cim, spec_metal);
			metalconf2statis_print(mcs);
		}
	}

/*
// AllPdb->output
	// from titles.txt, read every sdf's complex.pdb
	std::string pdbtitle = getp(parmap, "pdbs_titles");
	int num_th = std::stoi(getp(parmap, "num_th"));
	auto cims = getconfsfromtitlefile_pdb(pdbtitle, num_th);

	// new version.
	// statistic: bond_length, angle and dihedral angle.
	std::string spec_metal = getp(parmap, "specific_metal");
	if (spec_metal == "" || spec_metal == "ALL")
	{
		std::cout << "Statistics for all coordinating contacts." << std::endl;
		// extract metal and other atoms
		std::vector<MetalConf> mcs = aaconf2metalconf(cims);
		// release the memory
		vector<AAConformersInModel> ().swap(cims);
		metalconf2statis_print(mcs);
	}
	else
	{
		std::cout << "Statistics for " << spec_metal << std::endl;
		// extract metal and other atoms
		std::vector<MetalConf> mcs = aaconf2metalconf_spec(cims, spec_metal);
		// release the memory
		vector<AAConformersInModel> ().swap(cims);
		metalconf2statis_print(mcs);
	}
*/
}


