/*
 * scaffoldscreen.cpp
 *
 *  Created on: 2021年3月8日
 *      Author: yxchen
 */

/*
 * Screen Scaffolds, only keep those with e.g. monomor & <300aa.
 * Rank multiple scaffolds for the target scaffold (using their pocket volume).
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

int main(int argc, char **argv){
	string parfile(argv[1]);
	ifstream ifs;
	int MaxChainNum = 99; // maximum chain number in the target scaffold.
	int MaxAANum = 99999; // maximum amino acid number in each chain.
	int MinAANum = 0; // minium amino acid number in the longest chain.
	ScaPocket aimsp;
	vector<ScaPocket> poolsps;
	double radi = 0.0; // searching radius.
	double cubesize = 5.0; // side length of cube.
	string osd_path(""); // OnlyScreenDir
	bool NoMetal = false;
	bool NoApo = false;
	bool NoBig = false;
	int RankLigMW = 0; // ranking by ligand Molecule Weight (No hydrogen). (Simi -> diff)
	int InnerCA = 0; // if CA_num(<Radius) < InnerCA, then score + 100.

	ifs.open(parfile.c_str());
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
		if (words[0] == "MaxChainNum")
			MaxChainNum = stoi(words[2]);
		if (words[0] == "MaxAANum")
			MaxAANum = stoi(words[2]);
		if (words[0] == "MinAANum")
			MinAANum = stoi(words[2]);
		if (words[0] == "Radius")
			radi = stod(words[2]);
		if (words[0] == "CubeSize")
			cubesize = stod(words[2]);
		if (words[0] == "OnlyScreenDir" && words.size() == 3)
			osd_path = words[2];
		if (words[0] == "NoMetal" && words.size() == 3)
			if (words[2] == "Yes") NoMetal = true;
			else assert(words[2] == "No");
		if (words[0] == "NoApo" && words.size() == 3)
			if (words[2] == "Yes") NoApo = true;
			else assert(words[2] == "No");
		if (words[0] == "NoBig" && words.size() == 3)
			if (words[2] == "Yes") NoBig = true;
			else assert(words[2] == "No");
		if (words[0] == "RankLigMW" && words.size() == 3)
			RankLigMW = stoi(words[2]);
		if (words[0] == "InnerCA" && words.size() == 3)
			InnerCA = stoi(words[2]);
		if (words[0] == "AimPocket")
//			if (words.size() == 4)
//			{
//				string aimfile = words[2];
//				string ligname = words[3];
//				Pocket aimp;
//				aimp.readrefpdb(aimfile);
//				aimsp.poc2scapocket(aimfile, aimp, ligname, radi);
//			}
			if (words.size() == 6)
			{
				string aimfile = words[2];
				XYZ center {stod(words[3]), stod(words[4]), stod(words[5])};
				Pocket aimp;
				aimp.readrefpdb(aimfile);
				aimsp.poc2scapocket(aimfile, aimp, center, radi);
			}
			else
			{
				cout << "Wrong form of AimPocket." << endl;
				exit(1);
			}
		if (words[0] == "ScreenPocket")
//			if (words.size() == 4)
//			{
//				string scfile = words[2];
//				string ligname = words[3];
//				Pocket scp;
//				scp.readrefpdb(scfile);
//				if (passscreen(scp, MaxChainNum, MaxAANum, MinAANum, NoMetal, NoApo, NoBig))
//				{
//					ScaPocket sps;
//					sps.poc2scapocket(scfile, scp, ligname, radi);
//					poolsps.push_back(sps);
//				}
//			}
			if (words.size() == 6)
			{
				string scfile = words[2];
				XYZ center {stod(words[3]), stod(words[4]), stod(words[5])};
				Pocket scp;
				scp.readrefpdb(scfile);
				if (passscreen(scp, MaxChainNum, MaxAANum, MinAANum, NoMetal, NoApo, NoBig))
				{
					ScaPocket sps;
					sps.poc2scapocket(scfile, scp, center, radi);
					poolsps.push_back(sps);
				}
			}
			else
			{
				cout << "Wrong form of ScreenPocket." << endl;
				exit(1);
			}
	}
	ifs.close();

	vector<Pocket> poc_pools; // pocket after screen.
	vector<string> screen_pf; // pocket filename after screen.
	// screen by criteria.
	if (osd_path != "")
	{
		ofstream ofosd("screen.txt");
		cout << "Scaffolds screen begin:" << endl;
		vector<string> osd_files; // filenames
		getfile(osd_path, osd_files);
		for (auto pf : osd_files)
		{
			Pocket poc;
			poc.readrefpdb(pf);
			if (passscreen(poc, MaxChainNum, MaxAANum, MinAANum, NoMetal, NoApo, NoBig))
			{
				ofosd << pf << endl;
				poc_pools.push_back(poc);
				screen_pf.push_back(pf);
			}
		}
		ofosd.close();
		cout << "Screen End." << endl;
	}

	// 1. ranking by Ligand's MW
	if (RankLigMW > 0)
	{
		cout << "Ranking scaffolds by ligand's Molecule Weight: " << to_string(RankLigMW) << endl;
		map<int, vector<int>> rankmaps;
		if (InnerCA == 0)
			rankmaps = rankligmw(poc_pools, screen_pf, RankLigMW);
		else
			rankmaps = rankligmw(poc_pools, screen_pf, RankLigMW, InnerCA, radi);
		ofstream ofrlm("ranklmw.txt");
		for (auto iter = rankmaps.begin(); iter != rankmaps.end(); iter++)
			for (auto id : iter->second)
				ofrlm << iter->first << " " << screen_pf[id] << endl;
		cout << "Rank End." << endl;
//		exit(1);
	}
//	exit(1); // following hasn't been developed.

/*
	// calculate CubeNum
	int cn_aimsp = cubenumber(aimsp, cubesize);
	vector<int> cn_sps;
	for (auto psp : poolsps)
		cn_sps.push_back(cubenumber(psp, cubesize));

	// 2. ranked by similarity.
	vector<double> simis;
	map<double, vector<int>> simirank; // ranked by abs(simi), <poolsps_id>.
	int s = 0;
	for (auto cs : cn_sps)
	{
		double simi = (double)(cn_aimsp-cs)/cn_aimsp;
		simis.push_back(simi);
		simirank.insert(make_pair(fabs(simi), s));
		s++;
	}

	cout << "Ranking " << poolsps.size() << " scaffolds by pocket cavity similarity." << endl;
	for (auto iter = simirank.begin(); iter != simirank.end(); iter++)
		for (auto id : iter->second)
			cout << simis[id] << " " << poolsps[id].name() << endl;
*/
}
