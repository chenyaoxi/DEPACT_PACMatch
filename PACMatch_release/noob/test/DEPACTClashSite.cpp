/*
 * DEPACTClashSite.cpp
 *
 *  Created on: 2021年5月27日
 *      Author: yxchen
 */

#include "noob/theozyme.h"
#include "geometry/calculators.h"
#include <fstream>
#include <sstream>
using namespace NSPgeometry;
using namespace Theozyme;
using namespace std;

/*
 * Find clash sites for matched results.
 */

int main(int argc, char **argv)
{
	double pc_th = 2.7; // polar_clash_th.
	double npc_th = 3.2; // non-polar_clash_th.
	double ms_neigRMSD = 7.0; // RMSD threshold for finding ligand-Calpha neighbor site in MatchScaffold.
	double erotcutoff = 5.0;
	string filename = argv[1]; // input .params
	ifstream ifs;
	ifs.open(filename.c_str());
	if(!ifs.good()) {
		cout << ".params file failure" << endl;
		exit(1);
	}
	string mp_file, ms_file, ns_file; // MatchPocket, MatchScaffold, NativeScaffold.
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
		if (words[0] == "MatchPocket" && words.size() == 3)
			mp_file = words[2];
		if (words[0] == "MatchScaffold" && words.size() == 3)
			ms_file = words[2];
		if (words[0] == "MS_NeighborRMSD" && words.size() == 3)
			ms_neigRMSD = stod(words[2]);
		if (words[0] == "NativeScaffold" && words.size() == 3)
			ns_file = words[2];
		if (words[0] == "ERotCutoff" && words.size() == 3)
			erotcutoff = stod(words[2]);
	}
	if (mp_file == "" || ms_file == "")
	{
		cout << "input.par file should at least mentioned MatchPocket & MatchScaffold files." << endl;
		exit(1);
	}

	// output the input info.
	cout << "Input Info: " << endl;
	cout << "MatchPocket = " << mp_file << endl;
	cout << "MatchScaffold = " << ms_file << endl;
	cout << "MS_NeighborRMSD = " << ms_neigRMSD << endl;
	cout << "NativeScaffold = " << ns_file << endl;
	cout << endl;
	cout << "Results: " << endl;
	Pocket mp, ms, ns;
	mp.readrefpdb(mp_file);
	ms.readrefpdb(ms_file);
	if (ns_file != "")
		ns.readrefpdb(ns_file);
	auto mp_rs = mp.residues();
	auto mp_li = mp.ligands();
	auto ms_rs = ms.residues();

	// find NeighborSite in MatchScaffold for ligand
	cout << "Neighbor Sites of ligands (MatchPocket) in MatchScaffold are: " << endl;
	for (auto r : ms_rs)
		for (int i = 0; i < r.atmnames.size(); i++)
			if (r.atmnames[i] == "CA")
			{
				auto ca = r.atmcrds[i];
				bool neig = false;
				for (auto l : mp_li)
				{
					for (auto cl : l.atmcrds)
						if (ca.distance(cl) < ms_neigRMSD)
						{
							neig = true;
							break;
						}
					if (neig) break;
				}
				if (neig) cout << r.chainid << " " << r.resid << endl;
			}

	// find sites of MatchPocket in MatchScaffold.
	cout << "PocketResidues are in the sites of MatchScaffold: " << endl;
	map<int, int> mpsites; // mp_r,ms_site
	for (int i = 0; i < mp_rs.size(); i++)
		for (int j = 0; j < ms_rs.size(); j++)
		{
			auto pr = mp_rs[i];
			auto msr = ms_rs[j];
			int count = 0;
			for (auto prc : pr.atmcrds)
				for (auto msrc : msr.atmcrds)
					if (prc.x_ == msrc.x_ && prc.y_ == msrc.y_ && prc.z_ == msrc.z_)
						count++;
			if (count >= 4)
				mpsites.insert(make_pair(i,j));
		}
	for (auto iter = mpsites.begin(); iter != mpsites.end(); iter++)
	{
		auto pr = mp_rs[iter->first];
		auto msr = ms_rs[iter->second];
		cout << "PocketResidue " << pr.chainid << " " << pr.resid << " " << pr.resname << " matches in "
				<< msr.chainid << " " << msr.resid << endl;
	}

	// find clash sites between MatchPocket & MatchScaffold (considering SCRotamerLib).
	map<string, vector<int>> clashsites_ms = findclashRotLib(erotcutoff, ms_rs, mp_rs, mp_li, mpsites, pc_th, npc_th); // chainid, vector<rid>
	cout << "Clash sites in MatchScaffold:" << endl;
	if (clashsites_ms.size() == 0)
		cout << "No clash." << endl;
	else
	{
		for (auto iter = clashsites_ms.begin(); iter != clashsites_ms.end(); iter++)
		{
			cout << iter->first << " : ";
			for (auto s : iter->second)
				cout << s << " ";
			cout << endl;
		}
	}

	// if NewScaffold, find clash sites between MP & NS.
	if (ns_file != "")
	{
		auto ns_rs = ns.residues();
		map<string, vector<int>> clashsites_ns; // chainid, vector<rid>
		for (auto nsr : ns_rs)
		{
			bool clash = false;
 			for (int i = 0; i < mp_rs.size(); i++)
			{
				if (mpsites.count(i) == 0) continue; // non-matched residue don't calculate clash.
				clash = checkclash_ResSC(nsr, mp_rs[i], pc_th, npc_th);
				if (clash) break;
			}
			for (auto pl : mp_li)
			{
				for (int a = 0; a < nsr.atmcrds.size(); a++)
					for (int b = 0; b < pl.atmcrds.size(); b++)
					{
						double r = nsr.atmcrds[a].distance(pl.atmcrds[b]);
						if (r < pc_th)
							clash = true;
						else if (r < npc_th)
						{
							if (nsr.elementnames[a] == "C" && pl.elementnames[b] == "C") clash = true;
							if (nsr.elementnames[a] == "C" && pl.elementnames[b] == "S") clash = true;
							if (nsr.elementnames[a] == "S" && pl.elementnames[b] == "C") clash = true;
							if (nsr.elementnames[a] == "S" && pl.elementnames[b] == "S") clash = true;
						}
					}
				if (clash) break;
			}
			if (clash)
			{
				if (clashsites_ns.count(nsr.chainid) == 0)
					clashsites_ns.insert(make_pair(nsr.chainid, vector<int> {nsr.resid}));
				else
					clashsites_ns[nsr.chainid].push_back(nsr.resid);
			}
		}
		cout << "Clash sites in NativeScaffold:" << endl;
		if (clashsites_ns.size() == 0)
			cout << "No clash." << endl;
		else
		{
			for (auto iter = clashsites_ns.begin(); iter != clashsites_ns.end(); iter++)
			{
				cout << iter->first << " : ";
				for (auto s : iter->second)
					cout << s << " ";
				cout << endl;
			}
		}
	} // if consider nativescaffold.
}
