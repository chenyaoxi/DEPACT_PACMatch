/*
 * contactcode.h
 *
 *  Created on: 2019年9月29日
 *      Author: yxchen
 */

#ifndef NOOB_CONTACTCODE_H_
#define NOOB_CONTACTCODE_H_

#include "proteinrep/pdbrecord.h"
#include "proteinrep/pdbreader.h"
#include "geometry/xyz.h"
#include "noob/subsite.h"
#include <fstream>
#include <sstream>
#include <iterator>
#include <numeric>

using namespace NSPproteinrep;
using namespace NSPgeometry;

namespace ContactCode
{
typedef std::vector<std::vector<std::vector<PdbRecord>>> PdbRecordGroup;
typedef std::pair<std::string, std::string> ResAtm;
typedef std::pair<XYZ, std::vector<int>> ResCode;

struct SimCount
{
	int similar {0}; // atm num counted as similar
	int numa {0}; // atm num for A
	int numb {0}; // atm num for B
};

class LigAtmEnv
{
/* describe the env. for certain ligand atm. including ligatm's crd,
 * env's crd and env's type.
 */
public :
	void readligatmcrd(XYZ crd) {ligcrd_ = crd;}
	XYZ ligcrd() {return ligcrd_;}
	void readligatmid(int id) {ligatmid_ = id;}
	int ligatmid() {return ligatmid_;}
	void readligatmname(std::string name) {ligatmname_ = name;}
	std::string ligatmname() {return ligatmname_;}
	void buildrescodes(std::vector<ResCode> rcs)
	{
//		assert(rescodes_.size() == 0); // in some case we need to expand rescodes
		for (auto &rc : rcs)
			rescodes_.push_back(rc);
	}
	std::vector<ResCode> rescodes() {return rescodes_;}
//	void findctatmids(std::vector<int> ids)
//	{
//		assert(ctatmids_.size() == 0);
//		for (auto id : ids)
//			ctatmids_.push_back(id);
//	}
//	std::vector<int> ctatmids() {return ctatmids_;}
private :
	XYZ ligcrd_; // this lig-atm's crd
	int ligatmid_; // this lig-atm's id
	std::string ligatmname_; // this lig-atm's name
	std::vector<ResCode> rescodes_; // describe the env. of this lig-atm
//	std::vector<int> ctatmids_; // contain another 2 lig-atm ids; these three are in contacts
};

std::vector<LigAtmEnv> fetchpocenv_cb(PdbRecordGroup prg, std::map<ResAtm, std::vector<int>> contact_table);
SimCount calc_simcount(LigAtmEnv lae_a, LigAtmEnv lae_b, double rmsd_th);
std::vector<int> getcode(std::map<ResAtm, std::vector<int>> cmap, ResAtm ra);
std::map<std::vector<std::string>, std::vector<int>> readcodetable();
std::map<ResAtm, std::vector<int>> readcodetable_cavbase();
void readpocket(std::string prefix, int tot_num,
		std::vector<PdbRecordGroup> &poc_records);
std::string rname_3t1(std::string rname);
std::vector<int> encodecontact(std::string laname, std::string prname, std::string paname,
		std::map<std::vector<std::string>, std::vector<int>> contact_table);
void combinecode(std::vector<int> &outcode, std::vector<int> newcode);
std::vector<int> encodepocket(PdbRecordGroup p_r,
		std::map<std::vector<std::string>, std::vector<int>> contact_table);
LigAtmEnv fetchatmenv(PdbRecord la, PdbRecordGroup p_r, std::map<ResAtm, std::vector<int>> contact_table, int id);
LigAtmEnv fetchatmenv_cb(PdbRecord la, PdbRecordGroup p_r, std::map<ResAtm, std::vector<int>> contact_table, int id);
std::vector<LigAtmEnv> fetchpocenv(PdbRecordGroup p_r, std::map<ResAtm, std::vector<int>> contact_table);
std::vector<LigAtmEnv> fetchpocenv_partial(PdbRecordGroup p_r,
		std::map<ResAtm, std::vector<int>> contact_table, std::set<std::pair<std::string, std::string>> atmnames);
std::vector<int> encodepocket_partial(PdbRecordGroup p_r,
		std::map<std::vector<std::string>, std::vector<int>> contact_table, std::string par_clu);
bool incontact(std::string atmtype1, std::string atmtype2, double dist);
bool incontact_cb(PdbRecord la, PdbRecord pa, std::vector<PdbRecord> pr, std::map<ResAtm, std::vector<int>> contact_table);
int hamming_dist(std::vector<int> code_a, std::vector<int> code_b);
void hier_hamming(std::vector<std::vector<int>> poc_codes, double percent);
void fetchwholepocket(PdbRecordGroup ori_prg, PdbRecordGroup &new_prg, std::string ligname);
std::set<std::pair<std::string, std::string>> readpartialatm(std::string partial_cluster);
void output_simis_cavbase_all(std::vector<std::vector<LigAtmEnv>> poc_ligenvs, double rmsd_th);
void output_simis_cavbase_atom(std::vector<std::vector<LigAtmEnv>> poc_ligenvs, double rmsd_th, int poc_size);
std::vector<LigAtmEnv> combine_pool(std::vector<std::vector<LigAtmEnv>> poc_ligenvs);

struct AtmPair
{
	int i; // ligcrds_[ith]
	int j; // envcrds_[ith]
	std::vector<int> types;
//	bool operator < (const AtmPair &ap) const
//	{
//		if (i != ap.i) return i < ap.i;
//		else return j < ap.j;
//	}
};
class PocCode
{
public:
	void readcode(std::string filename);
	std::vector<XYZ> ligcrds() {return ligcrds_;}
	std::vector<XYZ> envcrds() {return envcrds_;}
	std::vector<AtmPair> contacts() {return contacts_;}
private:
	std::vector<XYZ> ligcrds_;
	std::vector<XYZ> envcrds_;
	std::vector<AtmPair> contacts_;
};
}


#endif /* NOOB_CONTACTCODE_H_ */
