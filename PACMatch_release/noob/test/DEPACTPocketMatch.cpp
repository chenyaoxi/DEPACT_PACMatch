/*
 * DEPACTPocketMatch.cpp
 *
 *  Created on: 2021年4月28日
 *      Author: yxchen
 */

#include "noob/theozyme.h"
#include "geometry/calculators.h"
#include <fstream>
#include <sstream>
#include <ctime>
#include <time.h>
using namespace NSPgeometry;
using namespace Theozyme;
using namespace std;
/*
 * Judge match for one-type pocket and one-type scaffolds (different confs.).
 */

int main(int argc, char **argv)
{
	clock_t startT, Tnow;
	startT = clock();
// read .params
	if (argc != 2)
	{
		cout << "The kam.par file is required as the only one input." << endl;
		exit(1);
	}
	string filename = argv[1]; // input .params
	string poc_path, sca_path, move_mode;
	string region = "";
	double RMSD_th = 2.0;
	string krfile = "";
	string atomlist = "";
	ifstream ifs;
	double plmc_clash = 1.5;
	double prpr_clash = 1.5;
	double ls_clash = 1.3;
	int keyresredconf = 1;
	int outputmatchnum = 100;
	double erotcutoff = 5.0;
	vector<int> rotlibc;
	ifs.open(filename.c_str());
	if(!ifs.good()) {
		cout << ".params file failure" << endl;
		exit(1);
	}
	while(true) {
		string line;
		line.clear();
		getline(ifs, line);
		if(!ifs.good()) break;
		if (line.size()==0) continue;
		vector<string> words;
		stringstream input(line);
		string word;
		while(input>>word) words.push_back(word);
		if (words[0] == "PocketDir" && words.size() == 3)
			poc_path = words[2];
		if (words[0] == "ScaffoldDir" && words.size() == 3)
			sca_path = words[2];
		if (words[0] == "MovePocket" && words.size() == 3)
			move_mode = words[2];
		if (words[0] == "Region" && words.size() == 3)
			region = words[2];
		if (words[0] == "RMSDCutoff" && words.size() == 3)
			RMSD_th = stod(words[2]);
		if (words[0] == "KeyResidues" && words.size() == 3)
			krfile = words[2];
		if (words[0] == "AtomList" && words.size() == 3)
			atomlist = words[2];
		if (words[0] == "LigScaMCClash" && words.size() == 3)
			plmc_clash = stod(words[2]);
		if (words[0] == "PocSCClash" && words.size() == 3)
			prpr_clash = stod(words[2]);
		if (words[0] == "LigPocSCClash" && words.size() == 3)
			ls_clash = stod(words[2]);
		if (words[0] == "RotLibChoice" && words.size() == 3)
			readrotlib(rotlibc, words[2]);
		if (words[0] == "KeyResRedConf" && words.size() == 3)
			keyresredconf = stoi(words[2]);
		if (words[0] == "OutputMatchNum" && words.size() == 3)
			outputmatchnum = stoi(words[2]);
		if (words[0] == "ERotCutoff" && words.size() == 3)
			erotcutoff = stod(words[2]);
	}
	ifs.close();

// ka_list
	AtomList al = readal(atomlist);
	Region reg = readregion(region);
// read pockets' ka_crd from pockets'
	// vector<string> filename -> vector<Pocket> -> vector<vector<KeyAtom>>
	vector<string> poc_files, sca_files; // filenames

	getfile(poc_path, poc_files);
	getfile(sca_path, sca_files);
	vector<Pocket> pockets;

	for (auto pf : poc_files)
	{
		Pocket poc;
		poc.readrefpdb(pf);
		pockets.push_back(poc);
	}

	vector<vector<KeyAtom>> poc_kas;
	for (auto &poc : pockets)
	{
		vector<KeyAtom> ka = extractka(al, poc);
		poc_kas.push_back(ka);
	}

	if (poc_kas.size() == 0)
	{
		cout << "PocketSize = 0, please check your PocketDir."  << endl;
		exit(1);
	}
	auto pks = poc_kas[0];
	if (pks.size() == 0)
	{
		cout << "Didn't get residues in the pocket." << endl;
		exit(1);
	}
	cout << "Pocket Residue Info:" << endl;
	for (int i = 0; i < pks.size(); i++)
		cout << " Pocket Residue " << i << ": " << pks[i].cid << " "
		    << pks[i].rid << " " << pks[i].rname << endl;

	vector<double> RMSD_ths(poc_kas[0].size(), RMSD_th);
	set<int> kr_ids;
	if (krfile != "")
		readkr_rmsdths(RMSD_ths, krfile, kr_ids);
	Tnow = clock();
	cout << "Time: " << (double)(Tnow-startT)/(3600*CLOCKS_PER_SEC) << " hours"
			<< " Finish reading pockets" << endl;
// for each scaffold, using ABACUS-like to predict possible vector<vector<KeyAtom>>
	int process = 0;
	std::map<string, vector<KeyAtom>> sca_kas; // filename_sca_ka
	std::map<double, MatchInfo> mis;
	std::map<string, vector<vector<XYZ>>> sca_crds_all; // filename_MCcrds.

	if (sca_files.size() == 0)
	{
		cout << "ScaffoldFiles = 0, please check your ScaffoldDir."  << endl;
		exit(1);
	}
	if (move_mode == "ON")
		cout << "Using GraphMatch." << endl;
	else if (move_mode == "OFF")
		cout << "Using DirectMatch." << endl;
	else
	{
		cout << "MovePocket should be ON or OFF." << endl;
		exit(1);
	}
	for (auto sf : sca_files)
	{
		process++;
		if (process % 100 == 0)
			cout << "have dealed with " << process << " scaffolds" << endl;
		Pocket sca;
		sca.readrefpdb(sf);
		readmcinfors_fromPocket(sf, sca, sca_crds_all);
		vector<KeyAtom> sca_ka;
		sca_ka = designka(erotcutoff, sca, al, reg, rotlibc);
		sca_kas.insert(make_pair(sf, sca_ka));
		std::cout << sf << " has " << sca_ka.size() << " possible rotamers." << std::endl;
		for (int pk = 0; pk < poc_kas.size(); pk++)
		{
			auto poc_ka = poc_kas[pk];
			if (move_mode == "OFF")
			{
				auto ligs = pockets[pk].ligands();
				// clash with mc?
				bool lig_ScaMC = false;
				for (auto scaCrds : sca_crds_all[sf])
				{
					for (auto scaCrd : scaCrds)
					{
						for (auto lig : ligs)
						{
							for (auto ligCrd : lig.atmcrds)
							{
								if (scaCrd.distance(ligCrd) < plmc_clash)
								{
									lig_ScaMC = true;
									break;
								}
								if (lig_ScaMC) break;
							}
							if (lig_ScaMC) break;
						}
						if (lig_ScaMC) break;
					}
					if (lig_ScaMC) break;
				}
				if (lig_ScaMC)
				{
					if (poc_files.size() == poc_kas.size())
						cout << "Ligands of Pocket " << poc_files[pk] << " have clash with ScaffoldMC of " << sf << endl;
					else
						cout << "Ligands of Pocket " << pk << " have clash with ScaffoldMC of " << sf << endl;
					continue;
				}
				directmatch_ka(poc_ka, sca_ka, RMSD_ths, mis, poc_files[pk], sf, kr_ids, keyresredconf);
			}
			else if (move_mode == "ON")
				graphicmatch_ka(poc_ka, sca_ka, RMSD_ths, mis, poc_files[pk], sf, kr_ids, keyresredconf);
		}
	}
	Tnow = clock();
	std::cout << "Time: " << (double)(Tnow-startT)/(3600*CLOCKS_PER_SEC) << " hours"
			<< " Finish matching with " << mis.size() << " matches." << std::endl;

//	// test
//	for (auto it = mis.begin(); it != mis.end(); it++)
//	{
//		auto mi = it->second;
//		cout << "Find: ";
//		for (auto v : mi.matchposes)
//			if (v != -1)
//				cout << sca_kas[mi.scaffold_filename][v].cid << " " << sca_kas[mi.scaffold_filename][v].rid << "; ";
//			else
//				cout << "Nan Nan; ";
//		cout << endl;
//	}

// combine two checking for speeding up
	std::map<double, MatchInfo> mis_new = double_check(outputmatchnum, mis, sca_crds_all, sca_kas, plmc_clash, prpr_clash, ls_clash, rotlibc);
	Tnow = clock();
	std::cout << "Time: " << (double)(Tnow-startT)/(3600*CLOCKS_PER_SEC) << " hours"
			<< " Finish checking and writing pockets with " << mis_new.size() << " non-redundant matches." << std::endl;

//	// test
//	for (auto it = mis_new.begin(); it != mis_new.end(); it++)
//	{
//		auto mi = it->second;
//		cout << "Find: ";
//		for (auto v : mi.matchposes)
//			if (v != -1)
//				cout << sca_kas[mi.scaffold_filename][v].cid << " " << sca_kas[mi.scaffold_filename][v].rid << "; ";
//			else
//				cout << "Nan Nan; ";
//		cout << endl;
//	}

// output mis
	ofstream ofs("match.txt");
	int i = 0;
	for (auto it = mis_new.begin(); it != mis_new.end(); it++)
	{
		if (i >= outputmatchnum) break;
		if (move_mode == "OFF")
			ofs << "Pocket " << it->second.pocket_filename
			<< " & Scaffold " << it->second.scaffold_filename
			<< " with score " << it->first << std::endl;
		else
		{
			std::string new_name = "new_Pocket_" + to_string(i) + ".pdb";
			ofs << "Pocket " << it->second.pocket_filename << " (" << new_name
			<< ") & Scaffold " << it->second.scaffold_filename
			<< " with score " << it->first << std::endl;
		}
		int num = 0;
		for (int pi = 0; pi < poc_kas[0].size(); pi++)
		{
			int si = it->second.matchposes[pi];
			string sf = it->second.scaffold_filename; // only being used to fetch the Res_info
			if (it->second.rmsds[pi] < 1.0) num++;
			if (si != -1)
				ofs << poc_kas[0][pi].cid << "_" << poc_kas[0][pi].rid << "_" << poc_kas[0][pi].rname << " : "
				<< sca_kas[sf][si].cid << "_" << sca_kas[sf][si].rid
				<< " with RMSD " << it->second.rmsds[pi] << std::endl;
			else
				ofs << poc_kas[0][pi].cid << "_" << poc_kas[0][pi].rid << "_" << poc_kas[0][pi].rname << " : "
				<< "no acceptable match" << std::endl;
		}
		if (num > 0)
			ofs << "Pocket " << it->second.pocket_filename
			<< " & Scaffold " << it->second.scaffold_filename
			<< " have good-match(<1.0A) " << std::to_string(num) << " in "
			<< std::to_string(poc_kas[0].size()) << std::endl;
		i++;
	}
	ofs.close();
	Tnow = clock();
	std::cout << "Time: " << (double)(Tnow-startT)/(3600*CLOCKS_PER_SEC) << " hours"
			<< " Finish the whole keyatommatch process." << std::endl;
}


