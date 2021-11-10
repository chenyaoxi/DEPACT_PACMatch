/*
 * readcrds.cpp
 *
 *  Created on: 2019年8月23日
 *      Author: yxchen
 */

// read crds and renew existed pdb file

#include "proteinrep/pdbreader.h"
#include "geometry/xyz.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
using namespace NSPproteinrep;

int main(int argc, char **argv)
{
	std::ifstream ifs(argv[1]);
	if (!ifs.good())
	{
		std::cout << "need crds file" << std::endl;
		exit(1);
	}
	std::vector<std::vector<double>> new_crds;
	while (true)
	{
		std::string line;
		line.clear();
		getline(ifs, line);
		if (!ifs.good()) break;
		if (line.size() == 0) continue;
		std::vector<double> ncrd;
		std::string word;
		std::stringstream input(line);
		while (input >> word) ncrd.push_back(std::stod(word));
		new_crds.push_back(ncrd);
	}


	int a = 0;
	PdbReader pr;
	pr.readpdb(argv[2]);
	std::string chainids = pr.chainids();
	std::ofstream ofs("new.pdb");
	for (int c = 0; c < chainids.size(); c++)
	{
		char cid = chainids[c];
		std::vector<std::string> seq = pr.getaminoacidsequence(cid);
		for (int r = 0; r < seq.size(); ++r)
		{
			typename PdbReader::ResKeyType reskey = pr.mappdbkeyint()->pdbResKey(r, c);
			std::vector<PdbRecord> &records = pr.records().at(cid).at(reskey);
			for (auto &record : records)
			{
				record.x = new_crds[a][0];
				record.y = new_crds[a][1];
				record.z = new_crds[a][2];
				a++;
				ofs << record.toString() << std::endl;
			}
		}
	}
}
