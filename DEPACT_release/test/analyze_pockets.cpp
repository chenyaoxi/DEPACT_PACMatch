/*
 * analyze_pockets.cpp
 *
 *  Created on: 2019年9月29日
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
	std::map<std::vector<std::string>, std::vector<int>> contact_table = readcodetable();

// read every pocket
	std::string prefix_pocket = getp(parmap, "Pocket");
	int num_poc(std::stoi(getp(parmap, "P_Num")));
	std::vector<PdbRecordGroup> poc_records;
	readpocket(prefix_pocket, num_poc, poc_records);

// partial cluster
	std::string par_clu = getp(parmap, "Partial_cluster");

// encode contact for every pocket
	std::vector<std::vector<int>> poc_codes; // [pocket_id][pocket_code]
	if (par_clu == "ALL")
		for (auto &p_r : poc_records)
			poc_codes.push_back(encodepocket(p_r, contact_table));
	else
		for (auto &p_r : poc_records)
			poc_codes.push_back(encodepocket_partial(p_r, contact_table, par_clu));
// read every native
	std::string prefix_native(getp(parmap, "Native"));
	int num_nat(std::stoi(getp(parmap, "N_Num")));
	std::vector<PdbRecordGroup> nat_records;
	readpocket(prefix_native, num_nat, nat_records);
// encode contact for every native
	if (par_clu == "ALL")
		for (auto &n_r : nat_records)
			poc_codes.push_back(encodepocket(n_r, contact_table));
	else
		for (auto &n_r : nat_records)
			poc_codes.push_back(encodepocket_partial(n_r, contact_table, par_clu));
//
//	// test
//	for (int i = 0; i < poc_codes.size(); i++)
//	{
//		std::cout << "ID_" << poc_ids[i] << " ";
//		for (auto pc : poc_codes[i])
//			std::cout << pc;
//		std::cout << std::endl;
//	}

// output codes
	std::string outcode = getp(parmap, "Output_code");
	if (outcode == "YES")
	{
		int id = 0;
		for (auto &pc : poc_codes)
		{
			std::cout << "Pocket_" << id++ << ": " << std::endl;
			for (auto & p : pc)
				std::cout << p;
			std::cout << std::endl;
		}
	}
	else if (outcode != "NO")
		std::cout << "Output_code should be YES or NO" << std::endl;

// cluster pockets by their codes
	double percent = std::stod(getp(parmap, "Percent"));
	hier_hamming(poc_codes, percent);
}
