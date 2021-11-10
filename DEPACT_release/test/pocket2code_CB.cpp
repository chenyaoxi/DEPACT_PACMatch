/*
 * pocket2code_CB.cpp
 *
 *  Created on: 2020年7月30日
 *      Author: yxchen
 */

#include "noob/contactcode.h"
#include "buildpocket.h"

using namespace NSPproteinrep;
using namespace ContactCode;
using namespace std;

/*
 * encode pocket by CB_criterion
 * input pocket.pdb
 * output pocket.pdb.txt
 */

int main(int argc, char **argv)
{
// read pocket.pdb.
	PdbReader pr;
	string poc(argv[1]);
	pr.readpdb(poc.c_str());
	std::string chainids = pr.chainids();
	PdbRecordGroup prg(chainids.size());
	for (int c = 0; c < chainids.size(); c++)
	{
		char cid = chainids[c];
		std::vector<std::string> seq = pr.getaminoacidsequence(cid);
		for (int r = 0; r < seq.size(); ++r)
		{
			typename PdbReader::ResKeyType reskey = pr.mappdbkeyint()->pdbResKey(r, c);
			std::vector<PdbRecord> &records = pr.records().at(cid).at(reskey);
			prg[c].push_back(records);
		}
	}

// read contactcodes.txt
	std::map<ResAtm, std::vector<int>> contact_table = readcodetable_cavbase();

// encode pocket.pdb. Default the first chain is for ligand.
	std::vector<LigAtmEnv> laes = fetchpocenv_cb(prg, contact_table);

// output codes. in pocket.pdb.txt
	string out = poc + ".txt";
	std::ofstream ofs(out.c_str());
	for (auto lae : laes)
	{
		auto ligcrd = lae.ligcrd();
		auto vrc = lae.rescodes();
		for (auto rc : vrc)
		{
			ofs << ligcrd.x_ << " " << ligcrd.y_ << " " << ligcrd.z_ << " "
					<< rc.first.x_ << " " << rc.first.y_ << " " << rc.first.z_;
			for (auto i : rc.second)
				ofs << " " << i;
			ofs << std::endl;
		} // for every ligand_atm
	}
	ofs.close();
}
