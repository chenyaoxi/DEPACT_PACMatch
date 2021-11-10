/*
 * contactcode.cpp
 *
 *  Created on: 2019年9月29日
 *      Author: yxchen
 */

#include "contactcode.h"
#include <dataio/datapaths.h>
using namespace NSPdataio;
using namespace NSPproteinrep;
using namespace ContactCode;
using namespace std;

std::map<std::vector<std::string>, std::vector<int>> ContactCode::readcodetable()
{
	std::map<std::vector<std::string>, std::vector<int>> codetable;
	std::ifstream ifcode("contactcodes.txt");
	if (!ifcode.good())
	{
		std::cout << "no contactcodes.txt" << std::endl;
		exit(1);
	}
	std::string line;
	getline(ifcode, line); // the first line is explanation.
	while(true)
	{
		line.clear();
		getline(ifcode, line);
		if (!ifcode.good()) break;
		if (line.size() == 0) continue;
		std::vector<std::string> lra; // ligand-residue-atom name
		std::vector<int> code; // e.g. 1 0 0 0 0 0 0 0
		std::vector<std::string> words;
		std::stringstream input(line);
		std::string word;
		while(input>>word) words.push_back(word);
		for (int i = 0; i < 3; i++)
			lra.push_back(words[i]);
		for (int i = 3; i < words.size(); i++)
			code.push_back(std::stoi(words[i]));
		codetable.insert(std::make_pair(lra, code));
	}
	ifcode.close();

	return codetable;
}

std::map<ResAtm, std::vector<int>> ContactCode::readcodetable_cavbase(){
	std::map<ResAtm, std::vector<int>> codetable;
	std::ifstream ifcode(getenvpath("PAAC_DATAPATH")+"contactcodes_CavBase.txt");
	if (!ifcode.good())
	{
		std::cout << "no contactcodes_CavBase.txt" << std::endl;
		exit(1);
	}
	std::string line;
	getline(ifcode, line); // the first line is explanation.
	while(true)
	{
		line.clear();
		getline(ifcode, line);
		if (!ifcode.good()) break;
		if (line.size() == 0) continue;
		ResAtm ra; // (resname, atmname)
		std::vector<int> code; // e.g. 4 or 4,5
		std::vector<std::string> words;
		std::stringstream input(line);
		std::string word;
		while(input>>word) words.push_back(word);
		ra.first = words[0];
		ra.second = words[1];
		for (int i = 2; i < words.size(); i++)
			code.push_back(std::stoi(words[i]));

		codetable.insert(std::make_pair(ra, code));
	}
	ifcode.close();

	return codetable;
}

std::string ContactCode::rname_3t1(std::string rname)
{
	if (rname == "GLY")  return "G";
	else if (rname == "ALA") return "A";
	else if (rname == "VAL") return "V";
	else if (rname == "LEU") return "L";
	else if (rname == "ILE") return "I";
	else if (rname == "MET") return "M";
	else if (rname == "TRP") return "W";
	else if (rname == "PHE") return "F";
	else if (rname == "PRO") return "P";
	else if (rname == "SER") return "S";
	else if (rname == "THR") return "T";
	else if (rname == "CYS") return "C";
	else if (rname == "TYR") return "Y";
	else if (rname == "ASN") return "N";
	else if (rname == "GLN") return "Q";
	else if (rname == "ASP") return "D";
	else if (rname == "GLU") return "E";
	else if (rname == "LYS") return "K";
	else if (rname == "ARG") return "R";
	else if (rname == "HIS") return "H";
	else if (rname == "HOH") return "X";
	else return "Z";
}

void ContactCode::readpocket(std::string prefix, int tot_num,
		std::vector<PdbRecordGroup> &poc_records)
{
//read every prefix_No.pdb in this dir. and return set{No., contacts}.
	for (int i = 0; i < tot_num; i++)
	{
		std::string fname = prefix+"_"+std::to_string(i)+".pdb";
		std::ifstream ifs(fname);
		if (!ifs.good()) continue;
		PdbReader pr;
		pr.readpdb(fname);
		std::string chainids = pr.chainids();
		std::vector<std::vector<std::vector<PdbRecord>>> pr_records(chainids.size());
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
		poc_records.push_back(pr_records);
	}
}

std::vector<int> ContactCode::encodecontact(std::string latype, std::string prname, std::string paname,
		std::map<std::vector<std::string>, std::vector<int>> contact_table)
{
	// laname: ligand_atm's type; prname: residue 3_name; paname: atm's full name.
	std::string pr1 = rname_3t1(prname);
	std::string patype = "X";
	patype[0] = paname[0];
	std::string xname = "X";
	std::vector<int> code(contact_table.begin()->second.size(), 0);
	std::vector<std::string> contact;
	// 1.0: residue_type independent
	contact.push_back(latype);
	contact.push_back(xname);
	contact.push_back(patype);
	if (contact_table.count(contact) > 0)
	{
		std::vector<int> newcode = contact_table[contact];
		combinecode(code, newcode);
	}
	contact.clear();
	// 1.5: residue_type independent & HOH(X)
	contact.push_back(latype);
	contact.push_back(xname);
	contact.push_back(xname);
	if (contact_table.count(contact) > 0)
	{
		std::vector<int> newcode = contact_table[contact];
		combinecode(code, newcode);
	}
	contact.clear();
	// 2.0: residue_type dependent & metal(Z)
	contact.push_back(xname);
	contact.push_back(pr1);
	contact.push_back(paname);
	if (contact_table.count(contact) > 0)
	{
		std::vector<int> newcode = contact_table[contact];
		combinecode(code, newcode);
	}
	contact.clear();

	return code;
}

void ContactCode::combinecode(std::vector<int> &outcode, std::vector<int> newcode)
{
	if (outcode.size() != newcode.size())
	{
		std::cout << "code length is different, cannot combine" << std::endl;
		exit(1);
	}
	for (int i = 0; i < outcode.size(); i++)
		outcode[i] = outcode[i] > newcode[i] ? outcode[i] : newcode[i];
}

std::vector<int> ContactCode::encodepocket(PdbRecordGroup p_r,
		std::map<std::vector<std::string>, std::vector<int>> contact_table)
{
	std::vector<int> poc_code;
	// default: first chain is ligand's chain.
	for (auto &lr : p_r[0])
		for (auto &la : lr)
		{
			std::vector<int> lacode(contact_table.begin()->second.size(), 0);
			for (int c = 1; c < p_r.size(); c++)
				for (auto &pr : p_r[c])
					for (auto &pa : pr)
					{
						std::string latype = la.namesymbol;
						std::string patype = pa.namesymbol;
						XYZ lacrd {la.x, la.y, la.z};
						double dist = lacrd.distance(XYZ {pa.x, pa.y, pa.z});
						if (!incontact(latype, patype, dist)) continue;
						std::vector<int> newcode = encodecontact(latype,
								pa.residuename, pa.atomname, contact_table);
						combinecode(lacode, newcode);
					} // every other atm
			poc_code.insert(poc_code.end(), lacode.begin(), lacode.end());
		} // every ligand atm
	return poc_code;
}

std::vector<int> ContactCode::getcode(std::map<ResAtm, std::vector<int>> cmap, ResAtm ra)
{
	std::vector<int> code;
	if (cmap.count(ra) == 1)
		code.assign(cmap[ra].begin(), cmap[ra].end());
	else
	{
		ResAtm pep(std::make_pair("PEP", ra.second));
		if (cmap.count(pep) == 1)
			code.assign(cmap[pep].begin(), cmap[pep].end());
	}
	if (code.size() == 0)
		code.push_back(-1);
	return code;
}

LigAtmEnv ContactCode::fetchatmenv(PdbRecord la, PdbRecordGroup p_r, std::map<ResAtm, std::vector<int>> contact_table, int id)
{
	LigAtmEnv lae;

	XYZ lacrd {la.x, la.y, la.z};
	lae.readligatmcrd(lacrd);
	lae.readligatmname(la.atomname);
	lae.readligatmid(id);

	std::vector<ResCode> rcs;
	for (int c = 1; c < p_r.size(); c++)
		for (auto &pr : p_r[c])
			for (auto &pa : pr)
			{
				std::string latype = la.namesymbol;
				std::string patype = pa.namesymbol;
				double dist = lacrd.distance(XYZ {pa.x, pa.y, pa.z});
				if (!incontact(latype, patype, dist)) continue;
				ResAtm ra(std::make_pair(pa.residuename, pa.atomname));
				std::vector<int> ra_code = getcode(contact_table, ra);
				if (ra_code[0] == -1) continue; // unknown atm
				ResCode rc(std::make_pair(XYZ {pa.x, pa.y, pa.z}, ra_code));
				rcs.push_back(rc);
			} // every poc atm
	if (rcs.size() > 0)
		lae.buildrescodes(rcs);

//	std::vector<int> ct_ids; // lig-atms. that contact with this lig-atm
//	int num = 0;
//	XYZ lbcrd;
//	for (auto &lr : p_r[0])
//	{
//		if (ct_ids.size() == 2) break;
//		for (auto &lb : lr)
//		{
//			double dist = lacrd.distance(XYZ {lb.x, lb.y, lb.z});
//			if (num != id && dist < 2.0)
//			{
//				ct_ids.push_back(num);
//				lbcrd.x_ = lb.x;
//				lbcrd.y_ = lb.y;
//				lbcrd.z_ = lb.z;
//			}
//			num++;
//			if (ct_ids.size() == 2) break;
//		}
//	}
//	if (ct_ids.size() == 1) // only 1 atm ct. with this lig-atm
//	{
//		id = ct_ids[0];
//		num = 0;
//		for (auto &lr : p_r[0])
//		{
//			if (ct_ids.size() == 2) break;
//			for (auto &lc : lr)
//			{
//				double dist = lbcrd.distance(XYZ {lc.x, lc.y, lc.z});
//				if (num != id && dist < 2.0)
//					ct_ids.push_back(num);
//				num++;
//				if (ct_ids.size() == 2) break;
//			}
//		}
//	}
//	lae.findctatmids(ct_ids);
	return lae;
}

std::vector<LigAtmEnv> ContactCode::fetchpocenv(PdbRecordGroup p_r, std::map<ResAtm, std::vector<int>> contact_table)
{
	std::vector<LigAtmEnv> laes;
	// default: first chain is ligand's chain.
	int id = 0;
	for (auto &lr : p_r[0])
		for (auto &la : lr)
			laes.push_back(fetchatmenv(la, p_r, contact_table, id++)); // every ligand atm
	return laes;
}

std::vector<LigAtmEnv> ContactCode::fetchpocenv_cb(PdbRecordGroup prg, std::map<ResAtm, std::vector<int>> contact_table)
{
	std::vector<LigAtmEnv> laes;
	// default: first chain is ligand's chain
	int id = 0;
	for (auto &lr : prg[0])
		for (auto &la : lr)
			laes.push_back(fetchatmenv_cb(la, prg, contact_table, id++));
	return laes;
}

LigAtmEnv ContactCode::fetchatmenv_cb(PdbRecord la, PdbRecordGroup p_r, std::map<ResAtm, std::vector<int>> contact_table, int id)
{
	LigAtmEnv lae;

	XYZ lacrd {la.x, la.y, la.z};
	lae.readligatmcrd(lacrd);
	lae.readligatmname(la.atomname);
	lae.readligatmid(id);

	std::vector<ResCode> rcs;
	for (int c = 1; c < p_r.size(); c++)
		for (auto &pr : p_r[c])
			for (auto &pa : pr)
			{
				std::string latype = la.namesymbol;
				std::string patype = pa.namesymbol;
//				double dist = lacrd.distance(XYZ {pa.x, pa.y, pa.z});
				if (!incontact_cb(la, pa, pr, contact_table)) continue;
				ResAtm ra(std::make_pair(pa.residuename, pa.atomname));
				std::vector<int> ra_code = getcode(contact_table, ra);
				if (ra_code[0] == -1) continue; // unknown atm
				ResCode rc(std::make_pair(XYZ {pa.x, pa.y, pa.z}, ra_code));
				rcs.push_back(rc);
			} // every poc atm
	if (rcs.size() > 0)
		lae.buildrescodes(rcs);
	return lae;
}

std::set<std::pair<std::string, std::string>> ContactCode::readpartialatm(std::string partial_cluster)
{
	// read partial_cluster.txt
	std::ifstream ifs(partial_cluster.c_str());
	if (!ifs.good())
	{
		std::cout << partial_cluster << " is not a file for atoms" << std::endl;
		exit(1);
	}
	std::set<std::pair<std::string, std::string>> atmnames;
	while(true)
	{
		std::string line;
		line.clear();
		getline(ifs, line);
		if (!ifs.good()) break;
		if (line.size() == 0) continue;
		if (line[0] == '#') continue;
		std::vector<std::string> words;
		std::stringstream input(line);
		std::string word;
		while (input >> word) words.push_back(word);
		atmnames.insert(std::make_pair(words[0], words[1]));
	}
	return atmnames;
}
std::vector<LigAtmEnv> ContactCode::fetchpocenv_partial(PdbRecordGroup p_r, std::map<ResAtm,
		std::vector<int>> contact_table, std::set<std::pair<std::string, std::string>> atmnames)
{


	std::vector<LigAtmEnv> laes;
	// default: first chain is ligand's chain.
	int id = 0;
	for (auto &lr : p_r[0])
		for (auto &la : lr)
			if (atmnames.count(std::make_pair(la.residuename, la.atomname)) == 1)
				laes.push_back(fetchatmenv(la, p_r, contact_table, id++)); // every ligand atm
	return laes;
}

std::vector<int> ContactCode::encodepocket_partial(PdbRecordGroup p_r,
		std::map<std::vector<std::string>, std::vector<int>> contact_table, std::string par_clu)
{
	// read partial group in ligand
	std::set<std::string> par_lig;
	std::ifstream ifs(par_clu.c_str());
	if (!ifs.good())
	{
		std::cout << par_clu << " is not good" << std::endl;
		exit(1);
	}
	while(true)
	{
		std::string line;
		line.clear();
		getline(ifs, line);
		if (!ifs.good()) break;
		if (line.size() == 1) continue;
		if (line[0] == '#') continue;
		par_lig.insert(line);
	}
	ifs.close();

	std::vector<int> poc_code;
	// default: first chain is ligand's chain.
	for (auto &lr : p_r[0])
		for (auto &la : lr)
		{
			if (par_lig.count(la.atomname) == 0) continue;
			std::vector<int> lacode(contact_table.begin()->second.size(), 0);
			for (int c = 1; c < p_r.size(); c++)
				for (auto &pr : p_r[c])
					for (auto &pa : pr)
					{
						std::string latype = la.namesymbol;
						std::string patype = pa.namesymbol;
						XYZ lacrd {la.x, la.y, la.z};
						double dist = lacrd.distance(XYZ {pa.x, pa.y, pa.z});
						if (!incontact(latype, patype, dist)) continue;
						std::vector<int> newcode = encodecontact(latype,
								pa.residuename, pa.atomname, contact_table);
						combinecode(lacode, newcode);
					} // every other atm
			poc_code.insert(poc_code.end(), lacode.begin(), lacode.end());
		} // every ligand atm
	return poc_code;
}

bool ContactCode::incontact(std::string atmtype1, std::string atmtype2, double dist)
{
	if (atmtype1 == "C" && atmtype2 == "C")
	{
		if (dist < 4.5) return true;
		else return false;
	}
	else if (atmtype2 == "S")
	{
		if (dist < 4.5) return true;
		else return false;
	}
	else
	{
		if (dist < 3.5) return true;
		else return false;
	}
}

bool ContactCode::incontact_cb(PdbRecord la, PdbRecord pa, std::vector<PdbRecord> pr, std::map<ResAtm, std::vector<int>> contact_table)
{
	// distance for all contacts
	std::string atmtype1 = la.namesymbol;
	std::string atmtype2 = pa.namesymbol;
	XYZ lacrd {la.x, la.y, la.z};
	XYZ pacrd {pa.x, pa.y, pa.z};
	double dist = lacrd.distance(pacrd);
	// non_residue: HOH & Metal
	if (rname_3t1(pa.residuename) == "X" || rname_3t1(pa.residuename) == "Z")
	{
		if (dist >= 3.0) return false;
		else return true;
	}
	// normal_residue
	else if (atmtype1 == "C" && atmtype2 == "C")
	{
		if (dist >= 4.5) return false;
		else return true;
	}
	else if (atmtype2 == "S")
	{
		if (dist >= 4.5) return false;
		else return true;
	}
	else
		if (dist >= 3.5) return false;

	// angle for the polar contact
	XYZ r {la.x-pa.x, la.y-pa.y, la.z-pa.z};
	// PEP {C-O} - first priority.
	if (pa.atomname == "O")
	{
		PdbRecord c;
		for (auto ar : pr)
			if (ar.atomname == "O")
				c = ar;
		XYZ v {pa.x-c.x, pa.y-c.y, pa.z-c.z};
		if (v*r/(v.length()*r.length()) > cos(100*M_PI/180)) return true;
		else return false;
	}
	// ARG {CD-NE-CZ} {CZ-NH1} {CZ-NH2}
	if (pa.residuename == "ARG")
	{
		PdbRecord cd, cz;
		for (auto ar : pr)
			if (ar.atomname == "CD")
				cd = ar;
			else if (ar.atomname == "CZ")
				cz = ar;
		if (pa.atomname == "NE")
		{
			XYZ v {2*pa.x-cd.x-cz.x, 2*pa.y-cd.y-cz.y, 2*pa.z-cd.z-cz.z};
			if (v*r/(v.length()*r.length()) > cos(100*M_PI/180)) return true;
			else return false;
		}
		else if (pa.atomname == "NH1" || pa.atomname == "NH2")
		{
			XYZ v {pa.x-cz.x, pa.y-cz.y, pa.z-cz.z};
			if (v*r/(v.length()*r.length()) > cos(100*M_PI/180)) return true;
			else return false;
		}
		else
			return true;
	}
	// ASN {CG-OD1} {CG-ND2}
	if (pa.residuename == "ASN")
	{
		PdbRecord cg;
		for (auto ar : pr)
			if (ar.atomname == "CG")
				cg = ar;
		if (pa.atomname == "OD1" || pa.atomname == "ND2")
		{
			XYZ v {pa.x-cg.x, pa.y-cg.y, pa.z-cg.z};
			if (v*r/(v.length()*r.length()) > cos(100*M_PI/180)) return true;
			else return false;
		}
		else
			return true;
	}
	// ASP {CG-OD1} {CG-OD2}
	if (pa.residuename == "ASP")
	{
		PdbRecord cg;
		for (auto ar : pr)
			if (ar.atomname == "CG")
				cg = ar;
		if (pa.atomname == "OD1" || pa.atomname == "OD2")
		{
			XYZ v {pa.x-cg.x, pa.y-cg.y, pa.z-cg.z};
			if (v*r/(v.length()*r.length()) > cos(100*M_PI/180)) return true;
			else return false;
		}
		else
			return true;
	}
	// GLN {CD-OE1} {CD-NE2}
	if (pa.residuename == "GLN")
	{
		PdbRecord cd;
		for (auto ar : pr)
			if (ar.atomname == "CD")
				cd = ar;
		if (pa.atomname == "OE1" || pa.atomname == "NE2")
		{
			XYZ v {pa.x-cd.x, pa.y-cd.y, pa.z-cd.z};
			if (v*r/(v.length()*r.length()) > cos(100*M_PI/180)) return true;
			else return false;
		}
		else
			return true;
	}
	// GLU {CD-OE1} {CD-OE2}
	if (pa.residuename == "GLU")
	{
		PdbRecord cd;
		for (auto ar : pr)
			if (ar.atomname == "CD")
				cd = ar;
		if (pa.atomname == "OE1" || pa.atomname == "OE2")
		{
			XYZ v {pa.x-cd.x, pa.y-cd.y, pa.z-cd.z};
			if (v*r/(v.length()*r.length()) > cos(100*M_PI/180)) return true;
			else return false;
		}
		else
			return true;
	}
	// HIS {CG-ND1-CE1} {CE1-NE2-CD2}
	if (pa.residuename == "HIS")
	{
		PdbRecord cg, ce1, cd2;
		for (auto ar : pr)
			if (ar.atomname == "CG")
				cg = ar;
			else if (ar.atomname == "CE1")
				ce1 = ar;
			else if (ar.atomname == "CD2")
				cd2 = ar;
		if (pa.atomname == "ND1")
		{
			XYZ v {2*pa.x-cg.x-ce1.x, 2*pa.y-cg.y-ce1.y, 2*pa.z-cg.z-ce1.z};
			if (v*r/(v.length()*r.length()) > cos(120*M_PI/180)) return true;
			else return false;
		}
		else if (pa.atomname == "NE2")
		{
			XYZ v {2*pa.x-cd2.x-ce1.x, 2*pa.y-cd2.y-ce1.y, 2*pa.z-cd2.z-ce1.z};
			if (v*r/(v.length()*r.length()) > cos(120*M_PI/180)) return true;
			else return false;
		}
		else
			return true;
	}
	// LYS {CE-NZ}
	if (pa.residuename == "LYS")
	{
		PdbRecord ce;
		for (auto ar : pr)
			if (ar.atomname == "CE")
				ce = ar;
		if (pa.atomname == "NZ")
		{
			XYZ v {pa.x-ce.x, pa.y-ce.y, pa.z-ce.z};
			if (v*r/(v.length()*r.length()) > cos(100*M_PI/180)) return true;
			else return false;
		}
		else
			return true;
	}
	// SER {CB-OG}
	if (pa.residuename == "SER")
	{
		PdbRecord cb;
		for (auto ar : pr)
			if (ar.atomname == "CB")
				cb = ar;
		if (pa.atomname == "OG")
		{
			XYZ v {pa.x-cb.x, pa.y-cb.y, pa.z-cb.z};
			if (v*r/(v.length()*r.length()) > cos(120*M_PI/180)) return true;
			else return false;
		}
		else
			return true;
	}
	// THR {CD1-NE1-CE2}
	if (pa.residuename == "THR")
	{
		PdbRecord cd1, ce2;
		for (auto ar : pr)
			if (ar.atomname == "CD1")
				cd1 = ar;
			else if (ar.atomname == "CE2")
				ce2 = ar;
		if (pa.atomname == "NE")
		{
			XYZ v {2*pa.x-cd1.x-ce2.x, 2*pa.y-cd1.y-ce2.y, 2*pa.z-cd1.z-ce2.z};
			if (v*r/(v.length()*r.length()) > cos(120*M_PI/180)) return true;
			else return false;
		}
		else
			return true;
	}
	// TYR {CZ-OH}
	if (pa.residuename == "TYR")
	{
		PdbRecord cz;
		for (auto ar : pr)
			if (ar.atomname == "CZ")
				cz = ar;
		if (pa.atomname == "OH")
		{
			XYZ v {pa.x-cz.x, pa.y-cz.y, pa.z-cz.z};
			if (v*r/(v.length()*r.length()) > cos(120*M_PI/180)) return true;
			else return false;
		}
		else
			return true;
	}

	return true;
}

int ContactCode::hamming_dist(std::vector<int> code_a, std::vector<int> code_b)
{
	int dist = 0;
	assert(code_a.size() == code_b.size());
	for (int i = 0; i < code_a.size(); i++)
		dist += abs(code_a[i]-code_b[i]);
	return dist;
}

void ContactCode::hier_hamming(std::vector<std::vector<int>> poc_codes, double percent)
{
	std::map<int, std::vector<int>> hm_map;
	for (int i = 0; i < poc_codes.size(); i++)
	{
		std::vector<int> member;
		member.push_back(i);
		hm_map.insert(std::make_pair(i, member));
	}

	int th = floor(poc_codes[0].size() * percent);
	int all_min_dist = 0;
	while(all_min_dist < th)
	{
		int min_dist = poc_codes[0].size();
		int min_a, min_b;
		for (auto it_a = hm_map.begin(); it_a != std::prev(hm_map.end()); it_a++)
		{
			for (auto it_b = std::next(it_a); it_b != hm_map.end(); it_b++)
			{
				// calc. max dist for it_a & it_b
				int max_dist = 0;
				for (auto & id_a : it_a->second)
					for (auto & id_b : it_b->second)
					{
						int dist = hamming_dist(poc_codes[id_a], poc_codes[id_b]);
						if (dist > max_dist)
							max_dist = dist;
					}
				// fresh the closest group
				if (max_dist < min_dist)
				{
					min_dist = max_dist;
					min_a = it_a->first;
					min_b = it_b->first;
				}
			}
		}
		if (min_dist >= all_min_dist)
			all_min_dist = min_dist;
		if (min_dist < th)
		{
			for (auto b : hm_map[min_b])
				hm_map[min_a].push_back(b);
			hm_map.erase(min_b);
		}
	}

	// test
	std::cout << "final min_dist: " << all_min_dist << std::endl;
	int id = 0;
	for (auto iter = hm_map.begin(); iter != hm_map.end(); iter++)
	{
		std::cout << id++ << ": ";
		for (auto mem : iter->second)
			std::cout << mem << " ";
		std::cout << std::endl;
	}
}

SimCount ContactCode::calc_simcount(LigAtmEnv lae_a, LigAtmEnv lae_b, double rmsd_th)
{
	SimCount sc;
	auto larc = lae_a.rescodes();
	auto lbrc = lae_b.rescodes();
	sc.numa = larc.size();
	sc.numb = lbrc.size();
	std::vector<int> matcha(sc.numa, 0);
	std::vector<int> matchb(sc.numb, 0);

	for (int a = 0; a < sc.numa; a++)
	{
		auto &rca = larc[a];
		for (int b = 0; b < sc.numb; b++)
		{
			auto &rcb = lbrc[b];
			if (rca.first.distance(rcb.first) < rmsd_th)
				for (auto codea : rca.second)
					for (auto codeb : rcb.second)
						if (codea == codeb)
						{
							matcha[a] = 1;
							matchb[b] = 1;
						}
		}
	}
	int suma = accumulate(matcha.begin(), matcha.end(), 0);
	int sumb = accumulate(matchb.begin(), matchb.end(), 0);
	sc.similar = suma < sumb ? suma : sumb;

	return sc;
}

void ContactCode::fetchwholepocket(PdbRecordGroup ori_prg, PdbRecordGroup &new_prg, std::string ligname)
{
	// contact only calc. lig-res & lig-media-res, no lig-coenzyme.
	int lc, lr;
	for (int c = 0; c < ori_prg.size(); c++)
		for (int r = 0; r < ori_prg[c].size(); r++)
			if (ori_prg[c][r][0].residuename == ligname)
			{
				lc = c;
				lr = r;
			}
	auto &lig_res = ori_prg[lc][lr];
	new_prg[0].push_back(lig_res);

	// direct contact
	std::vector<std::pair<int, int>> media; // to record mediators.
	for (int c = 0; c < ori_prg.size(); c++)
		for (int r = 0; r < ori_prg[c].size(); r++)
		{
			auto &ori_res = ori_prg[c][r];
			if (ori_res[0].residuename == ligname) continue;
			bool fetch = false;
			for (auto &pa : ori_res)
			{
				if (fetch) break;
				for (auto &la : lig_res)
				{
					XYZ crd_pa {pa.x, pa.y, pa.z};
					double dist = crd_pa.distance(XYZ {la.x, la.y, la.z});
					if (dist < 3.5 && ori_res.size() == 1) // a mediator
						media.push_back(std::make_pair(c, r));
					if (incontact(la.namesymbol, pa.namesymbol, dist))
					{
						fetch = true;
						break;
					}
				}
			}
			if (fetch) new_prg[1].push_back(ori_res);
		}

	std::vector<bool> ism(media.size(), false);
	// indirect contact (may keep atoms that are not mediators)
	for (int c = 0; c < ori_prg.size(); c++)
		for (int r = 0; r < ori_prg[c].size(); r++)
		{
			auto &ori_res = ori_prg[c][r];
			if (ori_res.size() == 1 || ori_res[0].residuename == ligname) continue;
			for (int m = 0; m < media.size(); m++)
			{
				auto &pair = media[m];
				auto &mediator = ori_prg[pair.first][pair.second];
				bool fetch = false;
				for (auto &pa : ori_res)
				{
					if (fetch) break;
					for (auto &ma : mediator)
					{
						XYZ crd_pa {pa.x, pa.y, pa.z};
						double dist = crd_pa.distance(XYZ {ma.x, ma.y, ma.z});
						if (incontact(ma.namesymbol, pa.namesymbol, dist))
						{
							fetch = true;
							break;
						}
					}
				}
				if (fetch)
				{
					new_prg[1].push_back(ori_res);
					ism[m] = true;
				}
			} // for every mediator
		}

	// dele. fake media
	for (int m = 0; m < ism.size(); m++)
		if (!ism[m])
			for (int n = 0; n < new_prg[1].size(); n++)
			{
				auto &np = new_prg[1][n];
				auto &op = ori_prg[media[m].first][media[m].second];
				if (op[0].chainid == np[0].chainid && op[0].residueid == np[0].residueid)
					auto iter = new_prg[1].erase(new_prg[1].begin()+n);
			}
}

void ContactCode::output_simis_cavbase_all(std::vector<std::vector<LigAtmEnv>> poc_ligenvs, double rmsd_th)
{
	std::map<std::pair<int, int>, std::vector<SimCount>> simicounts;
	for (int i = 0; i < poc_ligenvs.size()-1; i++)
		for (int j = i+1; j < poc_ligenvs.size(); j++)
		{
			std::pair<int, int> idpair {i, j};
			std::vector<SimCount> scounts;
			for (int k = 0; k < poc_ligenvs[i].size(); k++)
				scounts.push_back(calc_simcount(poc_ligenvs[i][k], poc_ligenvs[j][k], rmsd_th));
			simicounts.insert(std::make_pair(idpair, scounts));
		}

// output scores
	for (auto iter = simicounts.begin(); iter != simicounts.end(); iter++)
	{
		std::cout << iter->first.first << " " << iter->first.second << " ";
		double suma = 0.0;
		double sumb = 0.0;
		double sim = 0.0;
		for (auto &sc : iter->second)
		{
			suma += sc.numa;
			sumb += sc.numb;
			sim += sc.similar;
//			sum += (sc.similar + 0.001)/ (sc.numa + sc.numb - sc.similar+0.001);
		}
		std::cout << sim/(suma+sumb-sim) << std::endl;
	}
}

void ContactCode::output_simis_cavbase_atom(std::vector<std::vector<LigAtmEnv>> poc_ligenvs, double rmsd_th, int poc_size)
{
	std::map<std::pair<int, int>, std::vector<SimCount>> simicounts;
	for (int i = 0; i < poc_size; i++)
		for (int j = poc_size; j < poc_ligenvs.size(); j++)
		{
			std::pair<int, int> idpair {i, j};
			std::vector<SimCount> scounts;
			for (int k = 0; k < poc_ligenvs[i].size(); k++)
				scounts.push_back(calc_simcount(poc_ligenvs[i][k], poc_ligenvs[j][k], rmsd_th));
			simicounts.insert(std::make_pair(idpair, scounts));
		}
// output scores
	std::vector<std::pair<std::pair<int, int>, double>> max_simis;
	double suma = 0.0;
	double sumb = 0.0;
	double sim = 0.0;
	for (auto iter = simicounts.begin(); iter != simicounts.end(); iter++)
	{
		if (max_simis.size() < iter->second.size())
		{
			std::pair<int, int> compair {iter->first.first, iter->first.second};
			for (auto &sc : iter->second)
				if (sc.numa + sc.numb == 0) max_simis.push_back(std::make_pair(compair, 1.0));
				else
				{
					suma = sc.numa;
					sumb = sc.numb;
					sim = sc.similar;
					max_simis.push_back(std::make_pair(compair, sim/(suma+sumb-sim)));
				}
		}
		else
		{
			for (int i = 0; i < iter->second.size(); i++)
			{
				auto &sc = iter->second[i];
				suma = sc.numa;
				sumb = sc.numb;
				sim = sc.similar;
				double new_sim = sim/(suma+sumb-sim);
				if (new_sim > max_simis[i].second)
				{
					max_simis[i].first.first = iter->first.first;
					max_simis[i].first.second = iter->first.second;
					max_simis[i].second = new_sim;
				}
			}
		} // every atom
	} // every pair
	// output max_simis
	for (auto &ms : max_simis)
		std::cout << ms.first.first << " " << ms.first.second << " " << ms.second << std::endl;
}

std::vector<LigAtmEnv> ContactCode::combine_pool(std::vector<std::vector<LigAtmEnv>> poc_ligenvs)
{
	std::vector<LigAtmEnv> pool_ligenvs;
	auto &poc1 = poc_ligenvs[0];
	for (auto &lae : poc1) pool_ligenvs.push_back(lae);
	for (int i = 1; i < poc_ligenvs.size(); i++)
	{
		auto &poc = poc_ligenvs[i];
		for (int j = 0; j < pool_ligenvs.size(); j++)
			pool_ligenvs[j].buildrescodes(poc[j].rescodes());
	}
	return pool_ligenvs;
}

void PocCode::readcode(std::string filename)
{
	ifstream ifs(filename.c_str());
	if (!ifs.good())
	{
		cout << filename << " is wrong" << std::endl;
		exit(1);
	}
	while (true)
	{
		string line;
		line.clear();
		getline(ifs, line);
		if (!ifs.good()) break;
		if (line.size() == 0) continue;
		vector<std::string> words;
		stringstream input(line);
		string word;
		while (input>>word) words.push_back(word);
		XYZ lc {stod(words[0]), stod(words[1]), stod(words[2])};
		XYZ ec {stod(words[3]), stod(words[4]), stod(words[5])};
		vector<int> types;
		for (int i = 6; i < words.size(); i++)
			types.push_back(stoi(words[i]));
		int ith = -1;
		int eth = -1;
		for (int i = 0; i < ligcrds_.size(); i++)
		{
			auto l = ligcrds_[i];
			if (l.x_ == lc.x_ && l.y_ == lc.y_ && l.z_ == lc.z_)
			{
				ith = i;
				break;
			}
		}
		if (ith == -1)
		{
			ith = ligcrds_.size();
			ligcrds_.push_back(lc);
		}
		for (int i = 0; i < envcrds_.size(); i++)
		{
			auto e = envcrds_[i];
			if (e.x_ == ec.x_ && e.y_ == ec.y_ && e.z_ == ec.z_)
			{
				eth = i;
				break;
			}
		}
		if (eth == -1)
		{
			eth = envcrds_.size();
			envcrds_.push_back(ec);
		}
		AtmPair ap;
		ap.i = ith;
		ap.j = eth;
		ap.types = types;
		contacts_.push_back(ap);
	}
	ifs.close();
}



