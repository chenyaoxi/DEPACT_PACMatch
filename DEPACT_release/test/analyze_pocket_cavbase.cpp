/*
 * analyze_pocket_cavbase.cpp
 *
 *  Created on: 2019年11月11日
 *      Author: yxchen
 */

#include "noob/contactcode.h"
#include "buildpocket.h"

using namespace NSPproteinrep;
using namespace ContactCode;
using namespace buildpocket;

// cluster and rank all pocket.pdb in one dir.
int main(int argc, char **argv)
{
// read par
	std::string parname(argv[1]);
	std::map<std::string, std::string> parmap = readpar(parname);

// read contactcodes.txt
	std::map<ResAtm, std::vector<int>> contact_table = readcodetable_cavbase();

// read every pocket
	std::string prefix_pocket = getp(parmap, "Pocket");
	int num_poc(std::stoi(getp(parmap, "P_Num")));
	std::vector<PdbRecordGroup> poc_records;
	readpocket(prefix_pocket, num_poc, poc_records);

// encode env. for every ligand in every designed pocket.
	std::string partial_cluster = getp(parmap, "Partial_cluster");
	std::set<std::pair<std::string, std::string>> atmnames;
	if (partial_cluster != "ALL")
		atmnames = readpartialatm(partial_cluster);
	std::vector<std::vector<LigAtmEnv>> poc_ligenvs; // [pocket_id][lig_atm]'s env.
	for (auto &p_r : poc_records)
		if (partial_cluster == "ALL")
			poc_ligenvs.push_back(fetchpocenv(p_r, contact_table));
		else
			poc_ligenvs.push_back(fetchpocenv_partial(p_r, contact_table, atmnames));
	int poc_size = poc_ligenvs.size(); // to record how many Pockets are read.

// read every native
	std::string prefix_native(getp(parmap, "Native"));
	int num_nat(std::stoi(getp(parmap, "N_Num")));
	std::vector<PdbRecordGroup> nat_records;
	readpocket(prefix_native, num_nat, nat_records);

// encode env. for every ligand in every native pocket.
	for (auto &n_r : nat_records)
		if (partial_cluster == "ALL")
			poc_ligenvs.push_back(fetchpocenv(n_r, contact_table));
		else
			poc_ligenvs.push_back(fetchpocenv_partial(n_r, contact_table, atmnames));

// align and calc. similarity for every lig-atm.
	double rmsd_th(std::stod(getp(parmap, "RMSD_th")));
	std::string member(getp(parmap, "Compare_Member"));
	if (member == "ALL")
		output_simis_cavbase_all(poc_ligenvs, rmsd_th);
	else if (member == "ATOM")
		output_simis_cavbase_atom(poc_ligenvs, rmsd_th, poc_size);
	else
	{
		std::cout << "Compare_Member should be ALL or ATOM" << std::endl;
		exit(1);
	}
}
