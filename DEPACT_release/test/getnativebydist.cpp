/*
 * getnativebydist.cpp
 *
 *  Created on: 2019年10月8日
 *      Author: yxchen
 */

#include "noob/contactcode.h"
//#include "noob/subsite.h"

using namespace NSPproteinrep;
using namespace ContactCode;
//using namespace subsitedesign;

// get native pocket by contact_dist & calculate scores
int main(int argc, char **argv)
{
	// read native pdb
	std::string native_file(argv[1]);
	PdbReader pr;
	pr.readpdb(native_file);
	std::string chainids = pr.chainids();
	PdbRecordGroup pr_records(chainids.size());
	for (int c = 0; c < chainids.size(); c++)
	{
		char cid = chainids[c];
		std::vector<std::string> seq = pr.getaminoacidsequence(cid);
		for (int r = 0; r < seq.size(); ++r)
		{
			typename PdbReader::ResKeyType reskey = pr.mappdbkeyint()->pdbResKey(r, c);
			std::vector<PdbRecord> &records = pr.records().at(cid).at(reskey);
			pr_records[c].push_back(records);
		}
	}

	// fetch native pocket by contact_distance and output this pocket
	PdbRecordGroup wholepocket(2); // 0 for ligand & 1 for others
	std::string ligname(argv[2]);
	fetchwholepocket(pr_records, wholepocket, ligname);
	std::ofstream ofs("wholepocket.pdb");
	for (auto wpc : wholepocket)
		for (auto wpr : wpc)
			for (auto wpa : wpr)
				ofs << wpa.toString() << std::endl;
	ofs.close();

	// calculate scores
}
