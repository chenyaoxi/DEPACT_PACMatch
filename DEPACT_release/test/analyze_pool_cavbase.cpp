/*
 * analyze_pool_cavbase.cpp
 *
 *  Created on: 2019年11月26日
 *      Author: yxchen
 */

#include "noob/contactcode.h"
#include "buildpocket.h"

using namespace NSPproteinrep;
using namespace ContactCode;
using namespace buildpocket;

// to see if native_pool is covered by fssites_pool
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

// read every native
	std::string prefix_native(getp(parmap, "Native"));
	int num_nat(std::stoi(getp(parmap, "N_Num")));
	std::vector<PdbRecordGroup> nat_records;
	readpocket(prefix_native, num_nat, nat_records);

// encode env. for every ligand in every native pocket.
	std::vector<std::vector<LigAtmEnv>> nat_ligenvs; // [pocket_id][lig_atm]'s env.
	for (auto &n_r : nat_records)
		if (partial_cluster == "ALL")
			nat_ligenvs.push_back(fetchpocenv(n_r, contact_table));
		else
			nat_ligenvs.push_back(fetchpocenv_partial(n_r, contact_table, atmnames));

// combine each pool & calc. similarity
	double rmsd_th(std::stod(getp(parmap, "RMSD_th")));
	std::vector<LigAtmEnv> poc_pool = combine_pool(poc_ligenvs);
	std::vector<LigAtmEnv> nat_pool = combine_pool(nat_ligenvs);
	std::vector<SimCount> scounts;
	assert(poc_pool.size() == nat_pool.size());
	for (int i = 0; i < poc_pool.size(); i++)
		for (int j = 0; j < nat_pool.size(); j++)
			if (nat_pool[j].ligatmname() == poc_pool[i].ligatmname())
			{
				scounts.push_back(calc_simcount(poc_pool[i], nat_pool[j], rmsd_th));
				continue;
			}

// output simis for atm & all
	double scb = 0.0;
	double sim = 0.0;
	std::cout << "Every native atom's environment coverage (ranked by ids in Pocket ligand):" << std::endl;
	for (auto &sc : scounts)
	{
		double sb = sc.numb;
		double sm = sc.similar;
		scb += sb;
		sim += sm;
		if (sc.numb == 0) std::cout << 1.0 << " ";
		else std::cout << sm/sb << " ";
	}
	std::cout << std::endl;
	std::cout << "corresponding atom names: " << std::endl;
	for (int i = 0; i < poc_pool.size(); i++)
		std::cout << poc_pool[i].ligatmname() << " ";
	std::cout << std::endl;
	std::cout << "All Similarity: " << sim/scb << std::endl;
}


