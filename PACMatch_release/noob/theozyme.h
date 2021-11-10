/*
 * theozyme.h
 *
 *  Created on: 2020年6月19日
 *      Author: yxchen
 */

#ifndef NOOB_THEOZYME_H_
#define NOOB_THEOZYME_H_

/*
 * All infor about theozyme, only used for sampling from Dunbrack.database.
 */
#include <iostream>
#include <vector>
#include <map>
#include "geometry/xyz.h"
#include "designseq/Rotamer.h"
#include <algorithm>
//#include "cliquer/cliquer.h"
#include "geometry/quatfit.h"
#include <set>

using namespace NSPgeometry;
using namespace NSPdesignseq;
using namespace std;
namespace Theozyme
{

bool isresidue(std::string resname);

// It's the basic form about a general residue.
struct GeneralRes
{
	std::string resname;
	std::vector<std::string> atmnames;
	std::vector<XYZ> atmcrds;
	std::string chainid;
	int resid;
	std::vector<std::string> namesymbols;
	std::vector<std::string> namemodifiers;
	std::vector<std::string> elementnames;
};
XYZ getatmcrd(GeneralRes gr, std::string atmname);

// Describe GeneralRes move along a contact: Bond change
struct BondMove
{
	std::vector<std::string> GRname; // the name of the GeneralRes being rotated.
	std::vector<std::string> GRcid; // original chainID
	std::vector<int> GRrid; // original residueID
	std::vector<bool> GRisl;
	std::string atm1_cid;
	int atm1_rid;
	std::string atm1_name;
	bool isligand1 {false};
	std::string atm2_cid;
	int atm2_rid;
	std::string atm2_name;
	bool isligand2 {false};
	double bon_bin; // every move step
	int sample_num {2}; // add & delete for sample_num times separately.
};

// Describe GeneralRes rotates based on an axial: Dihedral rotation
struct RotateAxial
{
	std::vector<std::string> GRname; // the name of the GeneralRes being rotated.
	std::vector<std::string> GRcid; // original chainID
	std::vector<int> GRrid; // original residueID
	std::vector<bool> GRisl;
	std::string atm1_cid;
	int atm1_rid;
	std::string atm1_name;
	bool isligand1 {false};
	std::string atm2_cid;
	int atm2_rid;
	std::string atm2_name;
	bool isligand2 {false};
	double ang_begin {0};
	double ang_end {360};
	int sample_num {18}; // rotate every 20°
};

// Describe the interacting_angle change: Angle change
struct AngleTorsion
{
	std::vector<std::string> GRname; // the name of the GeneralRes being rotated.
	std::vector<std::string> GRcid; // original chainID
	std::vector<int> GRrid; // original residueID
	std::vector<bool> GRisl;
	std::string fa1_cid; // fix_atom1
	int fa1_rid;
	std::string fa1_name;
	bool isligand1 {false};
	std::string fa2_cid;
	int fa2_rid;
	std::string fa2_name;
	bool isligand2 {false};
	std::string ma_cid; // move_atom
	int ma_rid;
	std::string ma_name;
	bool isligand3 {false};
	double ang_bin; // every tor. change ang_bin.
	int sample_num {18}; // tor. num (*2 = positive + negative tor.)
};

struct BondJoint
{
	int ai; // from 1 ~ 3
	int aj; // from 4 ~ 6
	double b_b; // bond_length begin
	double b_e; // bond_length end
};

struct AngleJoint
{
	int ai; // from 1 ~ 6
	int aj; // from 1 ~ 6
	int ak; // from 4 ~ 6
	double a_b; // angle begin
	double a_e; // angle end
};

struct DihedJoint
{
	int ai; // from 1 ~ 6
	int aj; // from 1 ~ 6
	int ak; // from 1 ~ 6
	int al; // from 4 ~ 6
	double d_b; // dihedral begin
	double d_e; // dihedral end
};

// contains the relationship between Fixed_Parts(fa1, fa2, fa3) & Moved_Parts(ma1, ma2, ma3).
struct Joint
{
// Move parameters
	BondJoint bj;
	AngleJoint aj1;
	AngleJoint aj2;
	DihedJoint dj1;
	DihedJoint dj2;
	DihedJoint dj3;
// moved GeneralResidues
	std::vector<std::string> GRname; // the name of the GeneralRes being moved.
	std::vector<std::string> GRcid; // original chainID
	std::vector<int> GRrid; // original residueID
	std::vector<bool> GRisl; // is ligand?
// atoms for this joint in the <fa1, fa2, fa3, ma1, ma2, ma3> sequence.
	std::vector<std::string> anames;
	std::vector<std::string> acids;
	std::vector<int> arids;
	std::vector<bool> aisl;
	std::vector<XYZ> acrds;
};

// All infors & functions about theosite
class Pocket
{
public:
	void readrefpdb(const std::string &filename); // read ref_pdb
	void buildcenter(const std::string &filename); // build pocket by bc.gmy
	void readtheozyme(const std::string &filename); // read theozyme.thz
	void copyli(std::vector<GeneralRes> li) {for (auto a : li) ligands_.push_back(a);}
	void copyre(std::vector<GeneralRes> re) {for (auto a : re) residues_.push_back(a);}
	void crdchange(bool isligand, int gid, std::vector<XYZ> crdnew); // is ligand; which GR; new_crd
	std::vector<GeneralRes> ligands() const {return ligands_;}
	std::vector<GeneralRes> residues() const {return residues_;}
	void writepocket(std::string outfile);
	void transform(QuatFit qf);
private:
	std::vector<GeneralRes> ligands_;
	std::vector<GeneralRes> residues_;
};

struct ResidueSample
{
	std::string triname;
	std::string cid;
	int rid;
	std::vector<std::string> alignatms; // size >= 3.
};
std::vector<ResidueSample> init_rs(std::string filename);
void residuesample(std::vector<Pocket> &poc, ResidueSample rs);
void residuesample(std::vector<Pocket> &poc, ResidueSample rs, bool keepconf);
std::vector<Pocket> clashpocket(Pocket poc_origin, std::vector<Pocket> poc_new);
void clusterpocket(std::vector<Pocket> &poc_existing, std::vector<Pocket> poc_new);
std::vector<std::string> splitbyunderline(std::string line); // split A_B_C into vector{A,B,C};
XYZ zmatrix2XYZ(const XYZ rk, double b); // b = bond(rk-rl)
XYZ zmatrix2XYZ(const XYZ rk, double b, const XYZ rj, double theta); // -180<theta<180 = angle(lkj) as there is no phi
XYZ zmatrix2XYZ(const XYZ rk, double b, const XYZ rj, double theta, const XYZ ri, double phi);
// -180<phi<+180 = dihedral(lkji); 0<=theta<=180 = angle(lkj);
std::vector<Joint> readjoints(std::string filename);
void preparejointcrd(Pocket poc, Joint &joi);
void pocketchangebyjoint(Pocket &poc, Joint joi);
void pocketjointbond(Pocket &poc, Joint &j);
void pocketjointangle(Pocket &poc, Joint &j, int i);
void pocketjointdihedral(Pocket &poc, Joint &j, int i);
void samplepocketbyjoints(std::vector<Pocket> &pp, Pocket poc, std::vector<Joint> &jois, int seed, int samplenum);
void samplepocketbyjoints(std::vector<Pocket> &pp, Pocket poc, std::vector<Joint> &jois, int seed, int samplenum, bool keepconf);

struct ResMatch
{
	std::string poc_cid;
	int poc_rid;
	std::string poc_rname;
	std::vector<double> rmsds;
	std::vector<std::string> sca_cids;
	std::vector<int> sca_rids;
	std::vector<std::string> sca_rnames;
};

struct MatchInfo
{
	std::string pocket_filename;
	std::string scaffold_filename;
	std::vector<int> matchposes; // pocket_residue i find pos = matchposes[i] in scaffold(region).
	// if doesn't match ,matchpos = -1.
	std::vector<double> rmsds; // if doesn't match, rmsd = RMSD_th
	QuatFit qf; // for graphicmatch
	bool dm {false}; // if dm (directmatch), then qf is empty.
	set<string> scacodes; // <cid_rid> in scaffold, defined for checking not matching in one site.
};

struct AtomList
{
	vector<string> cids;
	vector<string> rids;
	vector<string> rname;
	vector<vector<string>> anames;
	vector<int> ids; // ids: distinguish same type res.
};

AtomList readal(string filename);

struct KeyAtom
{
	string rname;
	vector<string> anames;
	vector<XYZ> crds; // ka crds.
	int id; // distinguish same-type key residues.
	string cid;
	int rid;

	// the following elements are only used for Rotamer_rebuild.
	float phi;
	float psi;
	vector<XYZ> ncaco;
};

struct Region
{
	vector<string> cids;
	vector<int> i_begins;
	vector<int> i_ends;
};

void getfile(std::string path, std::vector<std::string> &filenames);
void readmcinfors(string filename, std::map<std::string, vector<vector<XYZ>>> &mc_crds, vector<string> &mc_crn);
void readmcinfors(string filename, std::map<std::string, vector<vector<XYZ>>> &mc_crds, vector<string> &mc_crn, std::string region_file);
void readmcinfors_fromPocket(string sf, Pocket sca, std::map<std::string, vector<vector<XYZ>>> &mc_crds);
void readallinfors(string filename, std::map<std::string, vector<vector<XYZ>>> &mc_crds, vector<string> &mc_crn,
		std::map<std::string, vector<vector<XYZ>>> &lig_crds, vector<string> &lig_crn);
bool hasedge(vector<XYZ> acrds, vector<XYZ> bcrds, double RMSD_th);
void graphs2mis(std::map<double, MatchInfo> &mis, vector<MatchInfo> &ms,
		std::string pfname, std::string sfname, vector<double> RMSD_ths, set<int> kr_ids);
std::map<double, MatchInfo> checkclash_mis(std::map<double, MatchInfo> mis,
		std::map<std::string, std::vector<std::vector<XYZ>>> sca_crds_all, vector<double> RMSD_ths);
std::map<double, MatchInfo> checkclash_lig_sc(std::map<double, MatchInfo> mis,
		std::map<string, vector<vector<XYZ>>> sca_all_crds);
std::map<double, MatchInfo> check_write_poc(std::map<double, MatchInfo> mis_n1,
		std::map<string, vector<KeyAtom>> sca_kas);
std::map<double, MatchInfo> double_check(int outputmatchnum, std::map<double, MatchInfo> mis,
		std::map<string, vector<vector<XYZ>>> sca_all_crds,
		std::map<string, vector<KeyAtom>> sca_kas,
		double plmc_clash, double prpr_clash, double ls_clash,
		vector<int> rotlibc = vector<int>());
vector<KeyAtom> designka(double erotcutoff, Pocket sca, AtomList al, Region reg, vector<int> rotlibc = vector<int>());
KeyAtom extractka_fromrot(AtomList al, Rotamer* rot);
vector<KeyAtom> extractka(AtomList al, Pocket poc);
void readkr_rmsdths(vector<double> &RMSD_ths, string filename);
void readkr_rmsdths(vector<double> &RMSD_ths, string filename, set<int> &kr_ids); // recording ids: key residues
char tri2one(string tri); // ResName in designseq.h has the same function
Region readregion(string region);
pair<float, float> getppfromgrs(GeneralRes pgr, GeneralRes gr, GeneralRes ngr);
vector<vector<XYZ>> get_mc_crds(vector<GeneralRes> residues);
bool mc_rg_hasclash(vector<XYZ> sc_crds, vector<vector<XYZ>> sca_mc_crds, int thissite);
void graphicmatch_ka(vector<KeyAtom> poc_ka, vector<KeyAtom> sca_ka,
		vector<double> RMSD_ths, std::map<double, MatchInfo> &mis, string pfilename, string sfilename,
		set<int> kr_ids, int keyresred); // newest version.
void directmatch_ka(vector<KeyAtom> poc_ka, vector<KeyAtom> sca_ka,
		vector<double> RMSD_ths, std::map<double, MatchInfo> &mis, string pfilename, string sfilename,
		set<int> kr_ids, int keyresredconf);
vector<vector<int>> directmatch_oka(vector<KeyAtom> poc_ka_new, vector<KeyAtom> sca_ka, vector<double> RMSD_ths,
		vector<int> p_s);
void bnb_ka(vector<KeyAtom> poc_ka, vector<KeyAtom> sca_ka,
		vector<double> RMSD_ths, std::map<double, MatchInfo> &mis, string pfilename, string sfilename);
void renew_residues(Pocket &soc, vector<KeyAtom> sca_ka, MatchInfo mi, vector<int> rotlibc = vector<int>());
void readrotlib(vector<int> &rotlibc, string filename);
void addmatchinfo(vector<MatchInfo> &ms, vector<KeyAtom> sca_ka, vector<vector<int>> pidsid, vector<vector<double>> ps_rmsd,
		int pid, set<int> kr_ids, vector<double> RMSD_ths);
void addmi_simple(vector<vector<int>> &p_s_grps, int p, vector<int> sites, vector<KeyAtom> sca_ka, int p_num);

// Calculate Pocket volume of pocket in scaffold.
class ScaPocket
{
public:
//	void poc2scapocket(string filename, Pocket poc, string ligname, double radi); // pocket -> center_ & poccas_.
	void poc2scapocket(string filename, Pocket poc, XYZ center, double radi); // by center XYZ
	string name() const {return name_;}
	XYZ center() const {return center_;}
	std::vector<XYZ> pocmcs() const {return pocmcs_;}
	double radi() const {return radi_;}
private:
	string name_; // pocketfile name.
	XYZ center_; // pocket (or ligand) center.
	std::vector<XYZ> pocmcs_; // surrounding MainChainAtm crds.
	double radi_; // searching radius;
};

// judge if scaffold can pass the screen criteria: ChainNum, MaxAANum & MinAANum.
bool passscreen(Pocket poc, int cn, int maxn, int minn);
// also may requiring no metal/having ligand/no DNA RNA membrane.
bool passscreen(Pocket poc, int cn, int maxn, int minn, bool NoMetal, bool NoApo, bool NoBig);
// ranking scaffolds by ligands' Molecule Weight.
map<int, vector<int>> rankligmw(vector<Pocket> pocs, vector<string> pfns, int lmw);
// if CA_num(<radius) < InnerCA, then score + 100.
map<int, vector<int>> rankligmw(vector<Pocket> pocs, vector<string> pfns, int lmw, int InnerCA, double Radius);
int cubenumber(ScaPocket sp, double cubesize);

// check if the closest atoms between a-b < d.
bool withindistance(GeneralRes a, GeneralRes b, double d);
// find clash between MatchPocket & MatchScaffold, considering RotLib
map<string, vector<int>> findclashRotLib(double erotcutoff, vector<GeneralRes> ms_rs, vector<GeneralRes> mp_rs,
		vector<GeneralRes> mp_li, map<int, int> mpsites, double pc_th, double npc_th);
// see if clash between two ResSC, using polar_clash_th & non-polar_clash_th
bool checkclash_ResSC(GeneralRes a, GeneralRes b, double pc_th, double npc_th);

}

#endif /* NOOB_THEOZYME_H_ */
