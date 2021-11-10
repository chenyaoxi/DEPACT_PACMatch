/*
 * pdb2seq.cpp
 *
 *  Created on: 2019年8月28日
 *      Author: yxchen
 */

#include "proteinrep/pdbreader.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>
using namespace NSPproteinrep;

int main(int argc, char **argv)
{
	std::map<std::string, std::string> t2o;
	t2o.insert(std::make_pair("GLY", "G"));
	t2o.insert(std::make_pair("ALA", "A"));
	t2o.insert(std::make_pair("VAL", "V"));
	t2o.insert(std::make_pair("LEU", "L"));
	t2o.insert(std::make_pair("ILE", "I"));
	t2o.insert(std::make_pair("MET", "M"));
	t2o.insert(std::make_pair("TRP", "W"));
	t2o.insert(std::make_pair("PHE", "F"));
	t2o.insert(std::make_pair("PRO", "P"));
	t2o.insert(std::make_pair("SER", "S"));
	t2o.insert(std::make_pair("THR", "T"));
	t2o.insert(std::make_pair("CYS", "C"));
	t2o.insert(std::make_pair("TYR", "Y"));
	t2o.insert(std::make_pair("ASN", "N"));
	t2o.insert(std::make_pair("GLN", "Q"));
	t2o.insert(std::make_pair("ASP", "D"));
	t2o.insert(std::make_pair("GLU", "E"));
	t2o.insert(std::make_pair("LYS", "K"));
	t2o.insert(std::make_pair("ARG", "R"));
	t2o.insert(std::make_pair("HIS", "H"));

	std::ofstream ofs(argv[2]);
	PdbReader pr;
	pr.readpdb(argv[1]);
	std::string chainids = pr.chainids();
	for (int c = 0; c < chainids.size(); c++)
	{
		char cid = chainids[c];
		std::vector<std::string> seq = pr.getaminoacidsequence(cid);
		for (int r = 0; r < seq.size(); ++r)
		{
			typename PdbReader::ResKeyType reskey = pr.mappdbkeyint()->pdbResKey(r, c);
			std::vector<PdbRecord> &records = pr.records().at(cid).at(reskey);
			ofs << t2o[records[0].residuename];
		}
		ofs << std::endl;
	}
}
