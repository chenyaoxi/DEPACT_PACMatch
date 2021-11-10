/*
 * theozyme.cpp
 *
 *  Created on: 2020年6月19日
 *      Author: yxchen
 */

/*
 * All infor about theozyme, only used for sampling from Dunbrack.database.
 */

#include "theozyme.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <string>
#include "geometry/rotation.h"
#include "proteinrep/pdbrecord.h"
#include "geometry/quatfit.h"
#include "designseq/RotamerLib.h"
#include "designseq/ProteinRep.h"
#include "geometry/calculators.h"
#include "dataio/inputlines.h"
//#include <boost/filesystem.hpp>
//#include <boost/algorithm/string.hpp>

#include "designseq/StructureInfo.h"
#include "designseq/DesignParameters.h"
#include "designseq/S1EnergyTable.h"
#include "designseq/S2MatrixFinder.h"
#include "designseq/DesignTemplate.h"

#include "noob/myclique.h"

#include <sys/types.h>
#include <dirent.h>
#include <numeric>

using namespace NSPgeometry;
using namespace NSPproteinrep;
using namespace NSPdesignseq;
using namespace Theozyme;
using namespace MyClique;
using namespace std;

static string EraseSpace(string s)
{
	const char ch = ' ';
	s.erase(s.find_last_not_of(" ")+1);
	s.erase(0, s.find_first_not_of(" "));
	return s;
}

// judge if is residue
bool Theozyme::isresidue(std::string resname)
{
	std::vector<std::string> residues {"GLY", "ALA", "VAL", "LEU",
	"ILE", "MET", "TRP", "PHE", "PRO", "SER", "THR", "CYS", "TYR",
	"ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS"};
	bool isr = false;
	if (std::find(residues.begin(), residues.end(), resname) != residues.end())
		isr = true;
	return isr;
}

// read theozyme structures from .pdb
void Pocket::readrefpdb(const std::string &filename)
{
	// new_new version
	const vector<int> fw { 6, 6, 2, 2, 1, 4, 1, 4, 1, 11, 8, 8, 6, 6, 10, 2 };
	vector<vector<string>> infors;
	NSPdataio::TextLines lines;
	lines.init(filename);
	for (auto & line : lines.lines())
	{
		if (line.substr(0, 6) == "ATOM  " || line.substr(0, 6) == "HETATM")
		{
			std::string eline =line+std::string(40,' ');
			vector<string> elines = NSPdataio::parseline(eline, fw);
			infors.push_back(elines);
		}
	}
	std::string chainid_pre = infors[0][6];
	int resid_pre = stoi(infors[0][7]);
	GeneralRes gr;
	gr.resname = infors[0][5];
	gr.resid = resid_pre;
	gr.chainid = chainid_pre;
	string namesymbol = EraseSpace(infors[0][2]); //boost::trim_copy(infors[0][2]);
	string namemodifier = EraseSpace(infors[0][3]); //boost::trim_copy(infors[0][3]);
	if (namesymbol[0] != 'H') // if (namesymbol[namesymbol.size()-1] != 'H')
	{
		gr.atmnames.push_back(namesymbol + namemodifier);
		gr.elementnames.push_back(infors[0][15]);
		gr.namesymbols.push_back(namesymbol);
		gr.namemodifiers.push_back(namemodifier);
		XYZ crd0 {stod(infors[0][9]), stod(infors[0][10]), stod(infors[0][11])};
		gr.atmcrds.push_back(crd0);
	}
	for (int i = 1; i < infors.size(); i++)
	{
		if (infors[i][6] != chainid_pre || stoi(infors[i][7]) != resid_pre)
		{
			chainid_pre = infors[i][6];
			resid_pre = stoi(infors[i][7]);
			if (isresidue(infors[i-1][5]))
				residues_.push_back(gr);
			else
				ligands_.push_back(gr);
			gr.atmcrds.clear();
			gr.atmnames.clear();
			gr.namesymbols.clear();
			gr.namemodifiers.clear();
			gr.elementnames.clear();
			gr.resname = infors[i][5];
			gr.chainid = chainid_pre;
			gr.resid = resid_pre;
		}
		namesymbol = EraseSpace(infors[i][2]); //boost::trim_copy(infors[i][2]);
		namemodifier = EraseSpace(infors[i][3]); //boost::trim_copy(infors[i][3]);
		if (namesymbol[0] != 'H') //if (namesymbol[namesymbol.size()-1] != 'H')
		{
			gr.atmnames.push_back(namesymbol + namemodifier);
			gr.elementnames.push_back(infors[i][15]);
			gr.namesymbols.push_back(namesymbol);
			gr.namemodifiers.push_back(namemodifier);
			XYZ crd {stod(infors[i][9]), stod(infors[i][10]), stod(infors[i][11])};
			gr.atmcrds.push_back(crd);
		}
	}

/* old version: may cause failure when numbers are too long.
	std::vector<std::vector<std::string>> infors;
	while(true) {
		std::string line;
		line.clear();
		getline(ifs, line);
		if(!ifs.good()) break;
		if (line.size()==0) continue;
		std::vector<std::string> words;
		stringstream input(line);
		std::string word;
		while(input>>word) words.push_back(word);
		if (words[0] == "HETATM" || words[0] == "ATOM")
			infors.push_back(words);
	}
	ifs.close();
	std::string chainid_pre = infors[0][4];
	int resid_pre = stoi(infors[0][5]);
	GeneralRes gr;
	gr.resname = infors[0][3];
	gr.atmnames.push_back(infors[0][2]);
	XYZ crd0 {stod(infors[0][6]), stod(infors[0][7]), stod(infors[0][8])};
	std::string ori_ele = infors[0].back(); // usually is 11
	int l_oe = ori_ele.size();
	std::string element;
	if (ori_ele[l_oe-1] == '+' || ori_ele[l_oe-1] == '-')
		element = ori_ele.substr(0, l_oe-2);
	else
		element = ori_ele;
	gr.elementnames.push_back(element);
	gr.namesymbols.push_back(element);
	gr.namemodifiers.push_back(infors[0][2].substr(element.size()));

	gr.atmcrds.push_back(crd0);
	gr.chainid = chainid_pre;
	gr.resid = resid_pre;
	for (int i = 1; i < infors.size(); i++)
	{
		if (infors[i][4] != chainid_pre || stoi(infors[i][5]) != resid_pre)
		{
			chainid_pre = infors[i][4];
			resid_pre = stoi(infors[i][5]);
			if (isresidue(infors[i-1][3]))
				residues_.push_back(gr);
			else
				ligands_.push_back(gr);
			gr.atmcrds.clear();
			gr.atmnames.clear();
			gr.namesymbols.clear();
			gr.namemodifiers.clear();
			gr.elementnames.clear();
			gr.resname = infors[i][3];
			gr.chainid = chainid_pre;
			gr.resid = resid_pre;
		}
		gr.atmnames.push_back(infors[i][2]);
		XYZ crd {stod(infors[i][6]), stod(infors[i][7]), stod(infors[i][8])};
		gr.atmcrds.push_back(crd);
		ori_ele = infors[i].back(); // usually is [11]
		l_oe = ori_ele.size();
		if (ori_ele[l_oe-1] == '+' || ori_ele[l_oe-1] == '-')
			element = ori_ele.substr(0, l_oe-2);
		else
			element = ori_ele;
		gr.elementnames.push_back(element);
		gr.namesymbols.push_back(element);
		gr.namemodifiers.push_back(infors[i][2].substr(element.size()));
	}
*/
	if (isresidue(gr.resname))
		residues_.push_back(gr);
	else
		ligands_.push_back(gr);
}

XYZ Theozyme::zmatrix2XYZ(const XYZ rk, double b) // b = bond(rk-rl)
{
	return rk + XYZ(b, 0.0, 0.0);
}
XYZ Theozyme::zmatrix2XYZ(const XYZ rk, double b, const XYZ rj, double theta)
{
	// -M_PI<theta<M_PI = angle(lkj), because there is no phi.
	XYZ rjk = rk - rj;
	try {
		XYZ ejk = rjk / std::sqrt(rjk.squarednorm());
		XYZ ei;
		double nxy = std::sqrt(ejk.x_ * ejk.x_ + ejk.y_ * ejk.y_);
		if (nxy > 0) {
			ei.x_ = ejk.y_ / nxy;
			ei.y_ = -ejk.x_ / nxy;
			ei.z_ = 0.0;
		} else {
			double nyz = std::sqrt(ejk.y_ * ejk.y_ + ejk.z_ * ejk.z_);
			ei.x_ = 0.0;
			ei.y_ = ejk.z_ / nyz;
			ei.z_ = -ejk.y_ / nyz;
		}
		XYZ rl = rk - (b * cos(theta*M_PI/180)) * ejk + (b * sin(theta*M_PI/180)) * ei;
		return rl;
	} catch (std::exception &e) {
		std::cout << "Error generating XYZ from zmatrix" << std::endl;
		exit(1);
	}
}
XYZ Theozyme::zmatrix2XYZ(const XYZ rk, double b, const XYZ rj, double theta, const XYZ ri, double phi)
{
	// 0<theta<M_PI = angle(lkj); -M_PI<phi<M_PI = dihedral(lkji)
	theta = abs(theta);
	if (theta > 180)
	{
		std::cout << theta << " is not in the region [0, 180]" << std::endl;
		exit(1);
	}
	XYZ rjk = rk - rj, rji = ri - rj;
	try {
		XYZ ejk = rjk / std::sqrt(rjk.squarednorm());
		rji= rji- dot(rji,ejk)*ejk;
		XYZ eji = rji/ std::sqrt(rji.squarednorm());
		XYZ rl_jk = - (b * cos(theta*M_PI/180)) * ejk;
		XYZ rl_ji = rk + (b * sin(theta*M_PI/180)) * eji; // > 0
		Rotation r;
		r.init(QuaternionCrd(ejk,phi),rk);
		r.apply(&rl_ji);
		XYZ rl = rl_jk + rl_ji ;
		return rl;
	} catch (std::exception &e) {
		std::cout << "Error generating XYZ from zmatrix" << std::endl;
		exit(1);
	}
}

std::vector<std::string> Theozyme::splitbyunderline(std::string line)
{
	std::vector<std::string> words;
	istringstream ist(line);
	string word;
	while (ist.good())
	{
		getline(ist, word, '_');
		if (!word.empty()) words.push_back(word);
	}
	return words;
}

void Pocket::buildcenter(const std::string &filename)
{
	// read file
	std::ifstream ifs;
	ifs.open(filename.c_str());
	if(!ifs.good()) {
		std::cout << ".gmy file failure" << std::endl;
		exit(1);
	}
	std::vector<std::vector<std::string>> infors;
	while(true) {
		std::string line;
		line.clear();
		getline(ifs, line);
		if(!ifs.good()) break;
		if (line.size()==0) continue;
		std::vector<std::string> words;
		stringstream input(line);
		std::string word;
		while(input>>word) words.push_back(word);
		infors.push_back(words);
	}
	ifs.close();

	// connect atms from residues_ and ligands
	std::map<std::string, XYZ> globalatms; // <cid_rid_atmname, crd>
	for (auto l : ligands_)
	{
		std::string name = l.chainid + "_" + to_string(l.resid) + "_";
		for (int i = 0; i < l.atmnames.size(); i++)
			globalatms.insert(std::make_pair(name+l.atmnames[i], l.atmcrds[i]));
	}
	for (auto r : residues_)
	{
		std::string name = r.chainid + "_" + to_string(r.resid) + "_";
		for (int i = 0; i < r.atmnames.size(); i++)
			globalatms.insert(std::make_pair(name+r.atmnames[i], r.atmcrds[i]));
	}

	// if init for three_atms
	bool init = false;
	if (infors[0][0] != "Residue" && infors[0][0] != "Ligand")
		init = true;
	if (init)
	{
		XYZ crdA{0.0, 0.0, 0.0};
		XYZ crdB = zmatrix2XYZ(crdA, stod(infors[1][2]));
		XYZ crdC;
		if (infors[2][1] == infors[0][0])
			crdC = zmatrix2XYZ(crdA, stod(infors[2][2]), crdB, stod(infors[2][4]));
		else
			crdC = zmatrix2XYZ(crdB, stod(infors[2][2]), crdA, stod(infors[2][4]));
		// has ref for ABC?
		if (infors[0][1] == "TO")
		{
			if (infors[1][3] != "TO" || infors[2][5] != "TO")
			{
				std::cout << "Error in three_atm_alignment" << std::endl;
				exit(1);
			}
			std::vector<std::string> aliA = splitbyunderline(infors[0][2]);
			std::vector<std::string> aliB = splitbyunderline(infors[1][4]);
			std::vector<std::string> aliC = splitbyunderline(infors[2][6]);
			std::vector<XYZ> crds_ori, crds_ali;
			XYZ crdaA, crdaB, crdaC;
			crds_ori.push_back(crdA); crds_ori.push_back(crdB); crds_ori.push_back(crdC);
			for (auto l : ligands_)
			{
				if (l.chainid == aliA[0] && to_string(l.resid) == aliA[1])
					crdaA = getatmcrd(l, aliA[2]);
				if (l.chainid == aliB[0] && to_string(l.resid) == aliB[1])
					crdaB = getatmcrd(l, aliB[2]);
				if (l.chainid == aliC[0] && to_string(l.resid) == aliC[1])
					crdaC = getatmcrd(l, aliC[2]);
			}
			for (auto r : residues_)
			{
				if (r.chainid == aliA[0] && to_string(r.resid) == aliA[1])
					crdaA = getatmcrd(r, aliA[2]);
				if (r.chainid == aliB[0] && to_string(r.resid) == aliB[1])
					crdaB = getatmcrd(r, aliB[2]);
				if (r.chainid == aliC[0] && to_string(r.resid) == aliC[1])
					crdaC = getatmcrd(r, aliC[2]);
			}
			crds_ali.push_back(crdaA); crds_ali.push_back(crdaB); crds_ali.push_back(crdaC);
			QuatFit qf;
			double rmsd = qf.fitting(crds_ali, crds_ori);
			crdA = crds_ori[0]; crdB = crds_ori[1]; crdC = crds_ori[2];
		} // do align
		globalatms.insert(std::make_pair(infors[0][0], crdA));
		globalatms.insert(std::make_pair(infors[1][0], crdB));
		globalatms.insert(std::make_pair(infors[2][0], crdC));
		std::cout << "Three axial atms' information:" << std::endl;
		std::cout << infors[0][0] << " " << crdA.x_ << " " << crdA.y_ << " " << crdA.z_ << std::endl;
		std::cout << infors[1][0] << " " << crdB.x_ << " " << crdB.y_ << " " << crdB.z_ << std::endl;
		std::cout << infors[2][0] << " " << crdC.x_ << " " << crdC.y_ << " " << crdC.z_ << std::endl;
	} // prepare for init

	// build center
	int start = 0;
	if (init) start = 3;
	for (int i = start; i < infors.size(); i++)
		if (infors[i][0] == "Residue")
		{
			GeneralRes gr;
			gr.chainid = infors[i][1];
			gr.resid = stoi(infors[i][2]);
			gr.resname = infors[i][3];
			// build rand rotamer
			XYZ ncrd {-45.590, 26.532, -32.562};
			XYZ cacrd {-45.080, 25.176, -32.806};
			XYZ ccrd {-45.934, 24.450, -33.835};
			XYZ ocrd {-47.134, 24.666, -33.902};
			RotamerLib* rotLib = new RotamerLib("A0");
			RotamerGroup* gp = rotLib->getAAGroup(gr.resname);
			Rotamer* rot = gp->rotList[0];
    		Residue* res = new Residue();
    		res->addAtom(new Atom("N", ncrd));
    		res->addAtom(new Atom("CA", cacrd));
    		res->addAtom(new Atom("C", ccrd));
    		res->addAtom(new Atom("O", ocrd));
    		res->buildRotamer(rot);
    		std::vector<std::string> resatms;
    		std::vector<XYZ> rescrds;
    		std::map<std::string, XYZ> atms_build; // <atom_name, crd>
    		for (int r = 0; r < res->getAtomList()->size(); r++)
    		{
    			auto atm = res->getAtomList()->at(r);
    			std::string fullname =  infors[i][1]+"_"+infors[i][2]+"_"+atm->name;
    			atms_build.insert(std::make_pair(fullname, atm->getCoord()));
    			resatms.push_back(atm->name);
    			rescrds.push_back(atm->getCoord());
    		}
			// infors -> align_atms
    		std::map<std::string, XYZ> aliatms;
			for (int j = 1; j <= stoi(infors[i][4]); j++)
			{
				XYZ cA = globalatms[infors[i+j][1]];
				double dist = stod(infors[i+j][2]);
				XYZ cB = globalatms[infors[i+j][3]];
				double ang = stod(infors[i+j][4]); // input is 0~180, but InternaltoXYZ need +-

				XYZ cC = globalatms[infors[i+j][5]];
				double dihed = stod(infors[i+j][6]);
				XYZ cn = zmatrix2XYZ(cA, dist, cB, ang, cC, dihed);
				std::string fullname = infors[i][1]+"_"+infors[i][2]+"_"+infors[i+j][0];
				globalatms.insert(std::make_pair(fullname, cn));
				aliatms.insert(std::make_pair(fullname, cn));
			}
			// do alignment and build gr
			std::vector<XYZ> alicrds, oricrds;
			for (auto it = aliatms.begin(); it != aliatms.end(); it++)
			{
				alicrds.push_back(it->second);
				oricrds.push_back(atms_build[it->first]);
			}
			QuatFit qf;
			qf.setup(alicrds, oricrds);
			qf.transform(rescrds);
			for (int j = 0; j < resatms.size(); j++)
			{
				gr.atmnames.push_back(resatms[j]);
				gr.atmcrds.push_back(rescrds[j]);
				std::string element;
				element.push_back(resatms[j][0]);
				gr.namesymbols.push_back(element);
				gr.elementnames.push_back(element);
				gr.namemodifiers.push_back(resatms[j].substr(element.size()));
				std::string fullname = infors[i][1]+"_"+infors[i][2]+"_"+resatms[j];
				if (globalatms.count(fullname) > 0)
					globalatms[fullname] = rescrds[j];
				else
					globalatms.insert(std::make_pair(fullname, rescrds[j]));
			}
			residues_.push_back(gr);
			i += stoi(infors[i][4]);
		}
		else if (infors[i][0] == "Ligand")
		{
			GeneralRes gr;
			gr.chainid = infors[i][1];
			gr.resid = stoi(infors[i][2]);
			gr.resname = infors[i][3];
			// build gr.
			for (int j = 1; j <= stoi(infors[i][4]); j++)
			{
				XYZ cA = globalatms[infors[i+j][1]];
				double dist = stod(infors[i+j][2]);
				XYZ cB = globalatms[infors[i+j][3]];
				double ang = stod(infors[i+j][4]);
				XYZ cC = globalatms[infors[i+j][5]];
				double dihed = stod(infors[i+j][6]); // to be 0~360
				XYZ cn = zmatrix2XYZ(cA, dist, cB, ang, cC, dihed);
				gr.atmnames.push_back(infors[i+j][0]);
				gr.atmcrds.push_back(cn);
				std::string element = infors[i+j][7];
				gr.namesymbols.push_back(element);
				gr.elementnames.push_back(element);
				gr.namemodifiers.push_back(infors[i+j][0].substr(element.size()));
				globalatms.insert(std::make_pair(infors[i][1]+"_"+infors[i][2]+"_"+infors[i+j][0], cn));
			}
			ligands_.push_back(gr);
			i += stoi(infors[i][4]);
		}
		else
		{
			std::cout << "Error in Residue/Ligand form" << std::endl;
			exit(1);
		}
}

void Pocket::readtheozyme(const std::string &filename)
{
	std::ifstream ifs;
	ifs.open(filename.c_str());
	if(!ifs.good()) {
		std::cout << ".thz file failure" << std::endl;
		exit(1);
	}
	std::string ref_pdb, bc_gmy; // ref_pdb & build center_geometry
	while(true) {
		std::string line;
		line.clear();
		getline(ifs, line);
		if(!ifs.good()) break;
		if (line.size()==0) continue;
		std::vector<std::string> words;
		stringstream input(line);
		std::string word;
		while(input>>word) words.push_back(word);
		if (words[0] == "ref_pdb" && words.size() == 3)
			ref_pdb = words[2];
		if (words[0] == "bc_geometry" && words.size() == 3)
			bc_gmy = words[2];
	}
	ifs.close();
	if (!ref_pdb.empty())
		readrefpdb(ref_pdb);
	if (!bc_gmy.empty())
		buildcenter(bc_gmy);
}

void Pocket::crdchange(bool isligand, int gid, std::vector<XYZ> crdnew)
{
	if (isligand)
	{
		ligands_[gid].atmcrds.clear();
		for (auto &c : crdnew)
			ligands_[gid].atmcrds.push_back(c);
	}
	else
	{
		residues_[gid].atmcrds.clear();
		for (auto &c : crdnew)
			residues_[gid].atmcrds.push_back(c);
	}
}

void Pocket::writepocket(std::string outfile)
{
	std::vector<PdbRecord> records;
	int aid = 1;
	int rid = 1;
	for (auto &p : ligands_)
	{
		for (int i = 0; i < p.atmnames.size(); i++)
		{
			PdbRecord record;
			record.label = "ATOM";
			record.atomname = p.atmnames[i];
			record.namesymbol = p.namesymbols[i];
			if (p.elementnames[i].size() == 1)
				record.elementname[1] = p.elementnames[i][0];
			else
			{
				record.elementname[0] = p.elementnames[i][0];
				record.elementname[1] = p.elementnames[i][1];
			}
			record.namemodifier = p.namemodifiers[i];
			record.residuename = p.resname;
			record.atomid = aid++;
			record.residueid = rid;
			record.x = p.atmcrds[i].x_;
			record.y = p.atmcrds[i].y_;
			record.z = p.atmcrds[i].z_;
			record.chainid = 'L';
			records.push_back(record);
		}
		rid++;
	}
	rid = 1;
	for (auto &p : residues_)
	{
		for (int i = 0; i < p.atmnames.size(); i++)
		{
			PdbRecord record;
			record.label = "ATOM";
			record.atomname = p.atmnames[i];
			record.namesymbol = p.namesymbols[i];
			if (p.elementnames[i].size() == 1)
				record.elementname[1] = p.elementnames[i][0];
			else
			{
				record.elementname[0] = p.elementnames[i][0];
				record.elementname[1] = p.elementnames[i][1];
			}
			record.namemodifier = p.namemodifiers[i];
			record.residuename = p.resname;
			record.atomid = aid++;
			record.residueid = rid;
			record.x = p.atmcrds[i].x_;
			record.y = p.atmcrds[i].y_;
			record.z = p.atmcrds[i].z_;
			record.chainid = 'R';
			records.push_back(record);
		}
		rid++;
	}

	std::ofstream ofs(outfile);
	for (auto &r : records) {
		ofs << r.toString()<<std::endl;
	}
	ofs.close();
}

void Pocket::transform(QuatFit qf)
{
	for (auto &l : ligands_)
		qf.transform(l.atmcrds);
	for (auto &r : residues_)
		qf.transform(r.atmcrds);
}

// initialize the std::vector<ResidueSample>
std::vector<ResidueSample> Theozyme::init_rs(std::string filename)
{
	std::vector<ResidueSample> rss;
	std::ifstream ifs(filename.c_str());
	if(!ifs.good()) {
		std::cout << ".smp file failure" << std::endl;
		exit(1);
	}
	std::vector<std::vector<std::string>> infors;
	while(true) {
		std::string line;
		line.clear();
		getline(ifs, line);
		if(!ifs.good()) break;
		if (line.size()==0) continue;
		std::vector<std::string> words;
		stringstream input(line);
		std::string word;
		while(input>>word) words.push_back(word);
		infors.push_back(words);
	}
	ifs.close();

	for (int i = 0; i < infors.size()/2; i++)
	{
		ResidueSample rs;
		rs.cid = infors[i*2][0];
		rs.rid = stoi(infors[i*2][1]);
		rs.triname = infors[i*2][2];
		for (auto atm : infors[i*2+1])
			rs.alignatms.push_back(atm);
		rss.push_back(rs);
	}
	return rss;
}

// read joints from file
std::vector<Joint> Theozyme::readjoints(std::string filename)
{
	std::ifstream ifs(filename.c_str());
	if(!ifs.good())
	{
		std::cout << ".rt file failure" << std::endl;
		exit(1);
	}
	std::vector<std::vector<std::string>> infors;
	while(true)
	{
		std::string line;
		line.clear();
		getline(ifs, line);
		if(!ifs.good()) break;
		if (line.size()==0) continue;
		std::vector<std::string> words;
		stringstream input(line);
		std::string word;
		while(input>>word) words.push_back(word);
		infors.push_back(words);
	}
	ifs.close();

	std::vector<Joint> js;
	for (int i = 0; i < infors.size(); i++)
	{

		if (infors[i][0] == "Joint:Begin")
		{
			if (infors.size() < i+11 || infors[i+10][0] != "Joint:End")
			{
				std::cout << "Joint in .rt is in wrong form." << std::endl;
				exit(1);
			}
			Joint j;
			assert(infors[i+1][0] == "MoveParts");
			for (int w = 0; w < (infors[i+1].size()-1)/4; w++)
			{
				j.GRcid.push_back(infors[i+1][1+w*4]);
				j.GRrid.push_back(stoi(infors[i+1][2+w*4]));
				j.GRname.push_back(infors[i+1][3+w*4]);
				if (infors[i+1][4+w*4] == "0")
					j.GRisl.push_back(false);
				else if (infors[i+1][4+w*4] == "1")
					j.GRisl.push_back(true);
				else
				{
					std::cout << "MoveParts: type should be 0 or 1" << std::endl;
					for (auto inf : infors[i+1])
						std::cout << inf << " ";
					std::cout << std::endl;
					exit(1);
				}
			}
			assert(infors[i+2][0] == "FixedAtoms");
			for (int w = 0; w < 3; w++)
			{
				j.acids.push_back(infors[i+2][1+w*4]);
				j.arids.push_back(stoi(infors[i+2][2+w*4]));
				j.anames.push_back(infors[i+2][3+w*4]);
				if (infors[i+2][4+w*4] == "0")
					j.aisl.push_back(false);
				else if (infors[i+2][4+w*4] == "1")
					j.aisl.push_back(true);
				else
				{
					std::cout << "FixedAtoms: type should be 0 or 1" << std::endl;
					for (auto inf : infors[i+2])
						std::cout << inf << " ";
					std::cout << std::endl;
					exit(1);
				}
			}
			assert(infors[i+3][0] == "MovedAtoms");
			for (int w = 0; w < 3; w++)
			{
				j.acids.push_back(infors[i+3][1+w*4]);
				j.arids.push_back(stoi(infors[i+3][2+w*4]));
				j.anames.push_back(infors[i+3][3+w*4]);
				if (infors[i+3][4+w*4] == "0")
					j.aisl.push_back(false);
				else if (infors[i+3][4+w*4] == "1")
					j.aisl.push_back(true);
				else
				{
					std::cout << "MovedAtoms: type should be 0 or 1" << std::endl;
					for (auto inf : infors[i+3])
						std::cout << inf << " ";
					std::cout << std::endl;
					exit(1);
				}
			}
			assert(infors[i+4][0] == "Bond");
			j.bj.ai = stoi(infors[i+4][1]);
			j.bj.aj = stoi(infors[i+4][2]);
			j.bj.b_b = stod(infors[i+4][3]);
			j.bj.b_e = stod(infors[i+4][4]);
			assert(infors[i+5][0] == "Angle");
			j.aj1.ai = stoi(infors[i+5][1]);
			j.aj1.aj = stoi(infors[i+5][2]);
			j.aj1.ak = stoi(infors[i+5][3]);
			j.aj1.a_b = stod(infors[i+5][4]);
			j.aj1.a_e = stod(infors[i+5][5]);
			assert(infors[i+6][0] == "Angle");
			j.aj2.ai = stoi(infors[i+6][1]);
			j.aj2.aj = stoi(infors[i+6][2]);
			j.aj2.ak = stoi(infors[i+6][3]);
			j.aj2.a_b = stod(infors[i+6][4]);
			j.aj2.a_e = stod(infors[i+6][5]);
			assert(infors[i+7][0] == "Dihedral");
			j.dj1.ai = stoi(infors[i+7][1]);
			j.dj1.aj = stoi(infors[i+7][2]);
			j.dj1.ak = stoi(infors[i+7][3]);
			j.dj1.al = stoi(infors[i+7][4]);
			j.dj1.d_b = stod(infors[i+7][5]);
			j.dj1.d_e = stod(infors[i+7][6]);
			assert(infors[i+8][0] == "Dihedral");
			j.dj2.ai = stoi(infors[i+8][1]);
			j.dj2.aj = stoi(infors[i+8][2]);
			j.dj2.ak = stoi(infors[i+8][3]);
			j.dj2.al = stoi(infors[i+8][4]);
			j.dj2.d_b = stod(infors[i+8][5]);
			j.dj2.d_e = stod(infors[i+8][6]);
			assert(infors[i+9][0] == "Dihedral");
			j.dj3.ai = stoi(infors[i+9][1]);
			j.dj3.aj = stoi(infors[i+9][2]);
			j.dj3.ak = stoi(infors[i+9][3]);
			j.dj3.al = stoi(infors[i+9][4]);
			j.dj3.d_b = stod(infors[i+9][5]);
			j.dj3.d_e = stod(infors[i+9][6]);
			js.push_back(j);
			i += 10;
		}
	}
	return js;
}

void Theozyme::preparejointcrd(Pocket poc, Joint &j)
{
	// prepare crdvalue from poc to joi.
	j.acrds.clear();
	for (int i = 0; i < 6; i++)
	{
		if (j.aisl[i])
		{
			auto ligands = poc.ligands();
			for (auto l : ligands)
				if (l.chainid == j.acids[i] && l.resid == j.arids[i])
					for (int a = 0; a < l.atmnames.size(); a++)
						if (l.atmnames[a] == j.anames[i])
							j.acrds.push_back(l.atmcrds[a]);
		}
		else
		{
			auto residues = poc.residues();
			for (auto r : residues)
				if (r.chainid == j.acids[i] && r.resid == j.arids[i])
					for (int a = 0; a < r.atmnames.size(); a++)
						if (r.atmnames[a] == j.anames[i])
							j.acrds.push_back(r.atmcrds[a]);

		}
	}
}

void Theozyme::pocketjointbond(Pocket &poc, Joint &j)
{
	preparejointcrd(poc, j);
	XYZ fa, ma;
	if (j.bj.ai <= 3)
	{
		fa = j.acrds[j.bj.ai-1];
		ma = j.acrds[j.bj.aj-1];
	}
	else
	{
		fa = j.acrds[j.bj.aj-1];
		ma = j.acrds[j.bj.ai-1];
	}
	double b_new;
	XYZ ma_move;
	if (j.bj.b_e != j.bj.b_b)
	{
		b_new = (j.bj.b_e-j.bj.b_b) * rand()/(double)(RAND_MAX) + j.bj.b_b;
		ma_move = (b_new/ma.distance(fa) - 1) * (ma - fa);
	}
	else
		ma_move = (j.bj.b_e - ma.distance(fa)) * (ma - fa);
	for (int g = 0; g < j.GRname.size(); g++)
	{
		std::vector<XYZ> crdnew;
		if (j.GRisl[g])
		{
			auto ligands = poc.ligands();
			for (int gid = 0; gid < ligands.size(); gid++)
				if (ligands[gid].chainid == j.GRcid[g] &&
						ligands[gid].resid == j.GRrid[g])
				{
					for (auto &co : ligands[gid].atmcrds)
						crdnew.push_back(co+ma_move);
					poc.crdchange(true, gid, crdnew);
					break;
				}
		}
		else
		{
			auto residues = poc.residues();
			for (int gid = 0; gid <residues.size(); gid++)
				if (residues[gid].chainid == j.GRcid[g] &&
						residues[gid].resid == j.GRrid[g])
				{
					for (auto &co : residues[gid].atmcrds)
						crdnew.push_back(co+ma_move);
					poc.crdchange(false, gid, crdnew);
					break;
				}
		}
	} // every GR
}

void Theozyme::pocketjointangle(Pocket &poc, Joint &j, int i)
{
	preparejointcrd(poc, j);
	AngleJoint aj;
	if (i == 1)
		aj = j.aj1;
	else
		aj = j.aj2;

// here fix and move is related to the angle_move here.
	XYZ fa1 = j.acrds[aj.ai-1];
	XYZ fa2 = j.acrds[aj.aj-1];
	XYZ ma = j.acrds[aj.ak-1];
	double a_change;
	XYZ axil = (ma-fa2) ^ (fa2-fa1);
	if (aj.a_e != aj.a_b)
	{
		double ang_new = (aj.a_e-aj.a_b) * rand()/(double)(RAND_MAX) + aj.a_b;
		a_change = ang_new - NSPgeometry::angle(fa1, fa2, ma)*180/M_PI;
	}
	else
		a_change = aj.a_e - NSPgeometry::angle(fa1, fa2, ma)*180/M_PI;
	Rotation rotat_pos(QuaternionCrd(axil, a_change), fa2);
	for (int g = 0; g < j.GRname.size(); g++)
	{
		std::vector<XYZ> crdnew;
		if (j.GRisl[g])
		{
			auto ligands = poc.ligands();
			for (int gid = 0; gid < ligands.size(); gid++)
				if (ligands[gid].chainid == j.GRcid[g] &&
						ligands[gid].resid == j.GRrid[g])
				{
					for (auto &co : ligands[gid].atmcrds)
						crdnew.push_back(rotat_pos.applytoCopy(co));
					poc.crdchange(true, gid, crdnew);
					break;
				}
		}
		else
		{
			auto residues = poc.residues();
			for (int gid = 0; gid <residues.size(); gid++)
				if (residues[gid].chainid == j.GRcid[g] &&
						residues[gid].resid == j.GRrid[g])
				{
					for (auto &co : residues[gid].atmcrds)
						crdnew.push_back(rotat_pos.applytoCopy(co));
					poc.crdchange(false, gid, crdnew);
					break;
				}
		}
	} // every generalresidue being rotated
}

void Theozyme::pocketjointdihedral(Pocket &poc, Joint &j, int i)
{
	preparejointcrd(poc, j);
	DihedJoint dj;
	if (i == 1)
		dj = j.dj1;
	else if (i == 2)
		dj = j.dj2;
	else
		dj = j.dj3;

// here fix and move is related to the dihedral_move here.
	XYZ fa1 = j.acrds[dj.aj-1];
	XYZ fa2 = j.acrds[dj.ak-1];
	XYZ ma = j.acrds[dj.al-1];
	double add, dihed_new;
	if (dj.d_e != dj.d_b)
	{
		if (dj.d_e > dj.d_b)
		    dihed_new = (dj.d_e - dj.d_b) * rand()/(double)(RAND_MAX) + dj.d_b;
		else
		{
			dihed_new = (dj.d_e + 180 + 180 - dj.d_b) * rand()/(double)(RAND_MAX) + dj.d_b;
			if (dihed_new > 180) dihed_new -= 360;
		}
	}
	else
		dihed_new = dj.d_e;
	double dihed_old = torsion(j.acrds[dj.ai-1], fa1, fa2, ma)*180/M_PI;
	if (dihed_old * dihed_new >= 0)
		add = dihed_new - dihed_old;
	else if (dihed_old < 0)
		add = - (dihed_old + 360 - dihed_new);
	else
		add = dihed_new - dihed_old + 360;
	Rotation rotat_pos(QuaternionCrd(fa2-fa1, add), fa2);
	for (int g = 0; g < j.GRname.size(); g++)
	{
		std::vector<XYZ> crdnew;
		if (j.GRisl[g])
		{
			auto ligands = poc.ligands();
			for (int gid = 0; gid < ligands.size(); gid++)
				if (ligands[gid].chainid == j.GRcid[g] &&
						ligands[gid].resid == j.GRrid[g])
				{
					for (auto &co : ligands[gid].atmcrds)
						crdnew.push_back(rotat_pos.applytoCopy(co));
					poc.crdchange(true, gid, crdnew);
					break;
				}
		}
		else
		{
			auto residues = poc.residues();
			for (int gid = 0; gid <residues.size(); gid++)
				if (residues[gid].chainid == j.GRcid[g] &&
						residues[gid].resid == j.GRrid[g])
				{
					for (auto &co : residues[gid].atmcrds)
						crdnew.push_back(rotat_pos.applytoCopy(co));
					poc.crdchange(false, gid, crdnew);
					break;
				}
		}
	} // every generalresidue being rotated
}

void Theozyme::pocketchangebyjoint(Pocket &poc, Joint j)
{
	// bond
	pocketjointbond(poc, j);
	// angle1
	pocketjointangle(poc, j, 1);
	// angle2
	pocketjointangle(poc, j, 2);
	// dihedral1
	pocketjointdihedral(poc, j ,1);
	// dihedral2
	pocketjointdihedral(poc, j ,2);
	// dihedral3
	pocketjointdihedral(poc, j ,3);
}

// sample pocket by joints
void Theozyme::samplepocketbyjoints(std::vector<Pocket> &pp, Pocket poc, std::vector<Joint> &jois, int seed, int samplenum)
{
	void srand(unsigned int seed);
	// do random sampling
	int no_change_num = 500;
	while (pp.size() < samplenum)
	{
		Pocket p_new = poc;
		for (auto j : jois)
			pocketchangebyjoint(p_new, j);
		std::vector<Pocket> p_temp {p_new};
		std::vector<Pocket> pp_noclash = clashpocket(poc, p_temp);
		if (pp_noclash.size() == 1)
		{
			pp.push_back(p_new);
			no_change_num = 500;
		}
		else
			no_change_num--;
		if (no_change_num == 0) break;
	}
}
void Theozyme::samplepocketbyjoints(std::vector<Pocket> &pp, Pocket poc, std::vector<Joint> &jois, int seed, int samplenum, bool keepconf)
{
	void srand(unsigned int seed);
	// do random sampling
	int no_change_num = 500;
	while (pp.size() < samplenum)
	{
		Pocket p_new = poc;
		for (auto j : jois)
			pocketchangebyjoint(p_new, j);
		if (keepconf)
		{
			std::vector<Pocket> p_temp {p_new};
			std::vector<Pocket> pp_noclash = clashpocket(poc, p_temp);
			if (pp_noclash.size() == 1)
			{
				pp.push_back(p_new);
				no_change_num = 500;
			}
			else
				no_change_num--;
		}
		else
			pp.push_back(p_new);
		if (no_change_num == 0) break;
	}
}

// do sampling from rot_Lib:A0 and B0 for every required residues
void Theozyme::residuesample(std::vector<Pocket> &poc, ResidueSample rs)
{
// fetch possible rotamer_confs.
// random choosed BBcrd
	XYZ ncrd {-45.590, 26.532, -32.562};
	XYZ cacrd {-45.080, 25.176, -32.806};
	XYZ ccrd {-45.934, 24.450, -33.835};
	XYZ ocrd {-47.134, 24.666, -33.902};
	vector<Rotamer*> rot_total; // mixed A1_lib and B1_lib
	RotamerLib* rotLib = new RotamerLib("A0");
	RotamerGroup* gp = rotLib->getAAGroup(rs.triname); // rotlists for certain res.
	for (auto rl : gp->rotList)
		rot_total.push_back(rl);
	rotLib = new RotamerLib("B0");
	gp = rotLib->getAAGroup(rs.triname);
	for (auto rl : gp->rotList)
		rot_total.push_back(rl);

// doing change
	std::vector<Pocket> new_poc;
	for (auto p : poc)
	{
		auto residues = p.residues();
		for (int gid = 0; gid < residues.size(); gid++)
		{
			auto r = residues[gid];
			if (r.chainid == rs.cid && r.resid == rs.rid)
			{
				std::vector<XYZ> crds_ali_old;
				for (auto an : rs.alignatms)
					for (int i = 0; i < r.atmnames.size(); i++)
						if (r.atmnames[i] == an)
							crds_ali_old.push_back(r.atmcrds[i]);
// for every possible rotamer_conf.
				for (auto rl : rot_total)
				{
		    		Residue* res = new Residue();
		    		res->addAtom(new Atom("N", ncrd));
		    		res->addAtom(new Atom("CA", cacrd));
		    		res->addAtom(new Atom("C", ccrd));
		    		res->addAtom(new Atom("O", ocrd));
		    		res->buildRotamer(rl);
		    		std::map<std::string, XYZ> atms_build; // <atom_name, crd>
		    		for (int i = 0; i < res->getAtomList()->size(); i++)
		    		{
		    			 auto atm = res->getAtomList()->at(i);
		    			atms_build.insert(std::make_pair(atm->name, atm->getCoord()));
		    		}
					std::vector<XYZ> crds_ali_new;
					for (auto an : rs.alignatms)
							crds_ali_new.push_back(atms_build[an]);
//align crds_ali_new -align2> crds_ali_old
					QuatFit qf;
					qf.setup(crds_ali_old, crds_ali_new);
// build crds_new from align.
					std::vector<XYZ> crds_new;
					for (auto a : r.atmnames)
						crds_new.push_back(atms_build[a]);
					qf.transform(crds_new);
//change crds_old -change2> crds_new
					Pocket p_smp;
					p_smp.copyli(p.ligands());
					p_smp.copyre(p.residues());
//					p_smp.copybm(p.bonmovs());
//					p_smp.copyro(p.rotaxls());
//					p_smp.copyan(p.angtors());
					p_smp.crdchange(false, gid, crds_new);
					std::vector<Pocket> poc_new {p_smp};
					std::vector<Pocket> poc_nc_new = clashpocket(poc[0], poc_new);
				//	clusterpocket(poc, poc_nc_new);
					for (auto pnn : poc_nc_new)
						new_poc.push_back(pnn);
				} // every rot_conf
				break; // this function only change one position's residue
			} // find the residue required conf_sample
		}
	} // every existing poc

	for (auto p : new_poc)
		poc.push_back(p);

}
void Theozyme::residuesample(std::vector<Pocket> &poc, ResidueSample rs, bool keepconf)
{
// fetch possible rotamer_confs.
// random choosed BBcrd
	XYZ ncrd {-45.590, 26.532, -32.562};
	XYZ cacrd {-45.080, 25.176, -32.806};
	XYZ ccrd {-45.934, 24.450, -33.835};
	XYZ ocrd {-47.134, 24.666, -33.902};
	vector<Rotamer*> rot_total; // mixed A0_lib and B0_lib
	RotamerLib* rotLib = new RotamerLib("A0");
	RotamerGroup* gp = rotLib->getAAGroup(rs.triname); // rotlists for certain res.
	for (auto rl : gp->rotList)
		rot_total.push_back(rl);
	rotLib = new RotamerLib("B0");
	gp = rotLib->getAAGroup(rs.triname);
	for (auto rl : gp->rotList)
		rot_total.push_back(rl);

// doing change & write
	std::vector<Pocket> new_poc;
	for (auto p : poc)
	{
		auto residues = p.residues();
		for (int gid = 0; gid < residues.size(); gid++)
		{
			auto r = residues[gid];
			if (r.chainid == rs.cid && r.resid == rs.rid)
			{
				std::vector<XYZ> crds_ali_old;
				for (auto an : rs.alignatms)
					for (int i = 0; i < r.atmnames.size(); i++)
						if (r.atmnames[i] == an)
							crds_ali_old.push_back(r.atmcrds[i]);
// for every possible rotamer_conf.
				for (auto rl : rot_total)
				{
		    		Residue* res = new Residue();
		    		res->addAtom(new Atom("N", ncrd));
		    		res->addAtom(new Atom("CA", cacrd));
		    		res->addAtom(new Atom("C", ccrd));
		    		res->addAtom(new Atom("O", ocrd));
		    		res->buildRotamer(rl);
		    		std::map<std::string, XYZ> atms_build; // <atom_name, crd>
		    		for (int i = 0; i < res->getAtomList()->size(); i++)
		    		{
		    			 auto atm = res->getAtomList()->at(i);
		    			atms_build.insert(std::make_pair(atm->name, atm->getCoord()));
		    		}
					std::vector<XYZ> crds_ali_new;
					for (auto an : rs.alignatms)
							crds_ali_new.push_back(atms_build[an]);
//align crds_ali_new -align2> crds_ali_old
					QuatFit qf;
					qf.setup(crds_ali_old, crds_ali_new);
// build crds_new from align.
					std::vector<XYZ> crds_new;
					for (auto a : r.atmnames)
						crds_new.push_back(atms_build[a]);
					qf.transform(crds_new);
//change crds_old -change2> crds_new
					Pocket p_smp;
					p_smp.copyli(p.ligands());
					p_smp.copyre(p.residues());
//					p_smp.copybm(p.bonmovs());
//					p_smp.copyro(p.rotaxls());
//					p_smp.copyan(p.angtors());
					p_smp.crdchange(false, gid, crds_new);
					vector<Pocket> psmp {p_smp};
					if (keepconf)
					{
						vector<Pocket> poc_nc_new = clashpocket(poc[0], psmp);
						for (auto pnn : poc_nc_new)
							new_poc.push_back(pnn);
					}
					else
					{
						new_poc.push_back(p_smp); // directly accepted.
//						vector<Pocket> poc_nc_new = clashpocket(p, psmp);
//						for (auto pnn : poc_nc_new)
//							new_poc.push_back(pnn);
					}
				} // every rot_conf
				break; // this function only change one position's residue
			} // find the residue required conf_sample
		}
	} // every existing poc

	for (auto p : new_poc)
		poc.push_back(p);

}

// get atm crd by atmname
XYZ Theozyme::getatmcrd(GeneralRes gr, std::string atmname)
{
	XYZ crd;
	bool find = false;
	for (int i = 0; i < gr.atmnames.size(); i++)
		if (gr.atmnames[i] == atmname)
		{
			crd.x_ = gr.atmcrds[i].x_;
			crd.y_ = gr.atmcrds[i].y_;
			crd.z_ = gr.atmcrds[i].z_;
			find = true;
		}
	if (!find)
	{
		std::cout << "cannot find " << atmname << " in " << gr.chainid << " "
				<< gr.resid << " " << gr.resname << std::endl;
		exit(1);
	}
	return crd;
}

// poc_new -clash filter> poc_new
std::vector<Pocket> Theozyme::clashpocket(Pocket poc_origin, std::vector<Pocket> poc_new)
{
	// simple_clash_filter: if poc_origin > 3.5A, then poc_new cannot < 3.5A
	std::map<int, std::vector<int>> with_th;
	double th = 3.5;
	int aid = 0;
	for (auto la : poc_origin.ligands())
		for (auto crda : la.atmcrds)
		{
			std::vector<int> neig_a;
			int bid = 0;
			for (auto lb : poc_origin.ligands())
				for (auto crdb : lb.atmcrds)
				{
					if (crda.distance(crdb) < 3.5)
						neig_a.push_back(bid);
					bid++;
				}
			for (auto rb : poc_origin.residues())
				for (auto crdb : rb.atmcrds)
				{
					if (crda.distance(crdb) < 3.5)
						neig_a.push_back(bid);
					bid++;
				}
			with_th.insert(make_pair(aid, neig_a));
			aid++;
		}
	for (auto ra : poc_origin.residues())
		for (auto crda : ra.atmcrds)
		{
			std::vector<int> neig_a;
			int bid = 0;
			for (auto lb : poc_origin.ligands())
				for (auto crdb : lb.atmcrds)
				{
					if (crda.distance(crdb) < 3.5)
						neig_a.push_back(bid);
					bid++;
				}
			for (auto rb : poc_origin.residues())
				for (auto crdb : rb.atmcrds)
				{
					if (crda.distance(crdb) < 3.5)
						neig_a.push_back(bid);
					bid++;
				}
			with_th.insert(make_pair(aid, neig_a));
			aid++;
		}

	// judge clash
	std::vector<Pocket> poc_nc; // no clash
	for (auto p : poc_new)
	{
		bool isnew = true;
		aid = 0;
		for (auto la : p.ligands())
			for (auto crda : la.atmcrds)
			{
				int bid = 0;
				for (auto lb : p.ligands())
					for (auto crdb : lb.atmcrds)
					{
						if (crda.distance(crdb) < 3.5)
						{
							auto neig_a = with_th[aid];
							auto it = find(neig_a.begin(), neig_a.end(), bid);
							if (it == neig_a.end())
								isnew = false;
						}
						bid++;
					}
				for (auto rb : p.residues())
					for (auto crdb : rb.atmcrds)
					{
						if (crda.distance(crdb) < 3.5)
						{
							auto neig_a = with_th[aid];
							auto it = find(neig_a.begin(), neig_a.end(), bid);
							if (it == neig_a.end())
								isnew = false;
						}
						bid++;
					}
				aid++;
			}
		for (auto ra : p.residues())
			for (auto crda : ra.atmcrds)
			{
				int bid = 0;
				for (auto lb : p.ligands())
					for (auto crdb : lb.atmcrds)
					{
						if (crda.distance(crdb) < 3.5)
						{
							auto neig_a = with_th[aid];
							auto it = find(neig_a.begin(), neig_a.end(), bid);
							if (it == neig_a.end())
								isnew = false;
						}
						bid++;
					}
				for (auto rb : p.residues())
					for (auto crdb : rb.atmcrds)
					{
						if (crda.distance(crdb) < 3.5)
						{
							auto neig_a = with_th[aid];
							auto it = find(neig_a.begin(), neig_a.end(), bid);
							if (it == neig_a.end())
								isnew = false;
						}
						bid++;
					}
				aid++;
			}
		if (isnew)
			poc_nc.push_back(p);
//		else
//			p.writepocket("clash_pocket.pdb");
	}
	return poc_nc;

}

// poc_new -cluster filter> poc_existing
void Theozyme::clusterpocket(std::vector<Pocket> &poc_existing, std::vector<Pocket> poc_new)
{
	// prepare MCcrds
	std::vector<std::vector<std::vector<XYZ>>> pn_mc_crds; // poc_res_mccrds
	std::vector<std::vector<std::vector<XYZ>>> po_mc_crds; // poc_res_mccrds
	for (auto po : poc_existing)
	{
		std::vector<std::vector<XYZ>> po_mcc;
		for (auto r : po.residues())
		{
			std::vector<XYZ> mc;
			mc.push_back(getatmcrd(r, "N"));
			mc.push_back(getatmcrd(r, "CA"));
			mc.push_back(getatmcrd(r, "C"));
			mc.push_back(getatmcrd(r, "O"));
			po_mcc.push_back(mc);
		}
		po_mc_crds.push_back(po_mcc);
	}
	for (auto pn : poc_new)
	{
		std::vector<std::vector<XYZ>> pn_mcc;
		for (auto r : pn.residues())
		{
			std::vector<XYZ> mc;
			mc.push_back(getatmcrd(r, "N"));
			mc.push_back(getatmcrd(r, "CA"));
			mc.push_back(getatmcrd(r, "C"));
			mc.push_back(getatmcrd(r, "O"));
			pn_mcc.push_back(mc);
		}
		pn_mc_crds.push_back(pn_mcc);
	}

	// cluster filter for a new pocket: at least one RMSD_res_mc > 1.0
	double RMSD_th = 1.0;
	std::vector<int> pass_ids;
	for (int i = 0; i < pn_mc_crds.size(); i++)
	{
		auto pnc  = pn_mc_crds[i];
		bool isnew = true;
		for (auto poc : po_mc_crds)
		{
			int count_small = 0;
			for (int r = 0; r < pnc.size(); r++)
				if (rmsd(pnc[r], poc[r]) < RMSD_th)
					count_small++;
			if (count_small == pnc.size())
				isnew = false;
		}
		if (isnew)
			pass_ids.push_back(i);
	}
	for (auto id : pass_ids)
		poc_existing.push_back(poc_new[id]);
}

// read all filenames from the dir's path
void Theozyme::getfile(std::string path, std::vector<std::string> &filenames)
{
    DIR *pDir;
    struct dirent* ptr;
    if(!(pDir = opendir(path.c_str())))
        return;
    while((ptr = readdir(pDir))!=0) {
        if (strcmp(ptr->d_name, ".") != 0 && strcmp(ptr->d_name, "..") != 0)
            filenames.push_back(path + "/" + ptr->d_name);
    }
    closedir(pDir);
}

// read all mc_crds from Pocket
void Theozyme::readmcinfors_fromPocket(string sf, Pocket sca, std::map<std::string, vector<vector<XYZ>>> &mc_crds)
{
	auto residues = sca.residues();
	vector<vector<XYZ>> ncacos;
	for (auto r : residues)
	{
		vector<XYZ> ncaco;
		ncaco.push_back(getatmcrd(r, "N"));
		ncaco.push_back(getatmcrd(r, "CA"));
		ncaco.push_back(getatmcrd(r, "C"));
		ncaco.push_back(getatmcrd(r, "O"));
		ncacos.push_back(ncaco);
	}
	mc_crds.insert(std::make_pair(sf, ncacos));
}


// read mainchain crds and related filenames. (for pocket)
void Theozyme::readmcinfors(string filename, std::map<std::string, vector<vector<XYZ>>> &mc_crds,
		vector<string> &mc_crn)
{
// fetch all mc
	Pocket poc;
	poc.readrefpdb(filename);

	auto residues = poc.residues();
	std::vector<std::vector<XYZ>> ncacs;
	for (int i = 0; i < residues.size(); i++)
	{
		auto r = residues[i];
		std::vector<XYZ> ncac;
		ncac.push_back(getatmcrd(r, "N"));
		ncac.push_back(getatmcrd(r, "CA"));
		ncac.push_back(getatmcrd(r, "C"));
//		ncac.push_back(getatmcrd(r, "O")); // as Ocrds is arbitrary for psi
		ncacs.push_back(ncac);
		if (mc_crn.size() < residues.size())
			mc_crn.push_back(r.chainid+"_"+to_string(r.resid)+"_"+r.resname);
		else
			assert(mc_crn[i] == r.chainid+"_"+to_string(r.resid)+"_"+r.resname);
	}
	mc_crds.insert(std::make_pair(filename, ncacs));
}
// read mainchain crds and related filenames. With regions (for scaffold)
void Theozyme::readmcinfors(string filename, std::map<std::string, vector<vector<XYZ>>> &mc_crds,
		vector<string> &mc_crn, std::string region_file)
{
// fetch region
	std::map<std::string, std::vector<int>> regions; // cid, {rid1, rid2, ...}
	ifstream ifs(region_file.c_str());
	if(!ifs.good()) {
		std::cout << ".reg file failure" << std::endl;
		exit(1);
	}
	while(true) {
		std::string line;
		line.clear();
		getline(ifs, line);
		if(!ifs.good()) break;
		if (line.size()==0) continue;
		std::vector<std::string> words;
		std::stringstream input(line);
		std::string word;
		while(input>>word) words.push_back(word);
		if (regions.count(words[0]) == 0)
		{
			std::vector<int> reg;
			for (int i = stoi(words[1]); i <= stoi(words[2]); i++)
				reg.push_back(i);
			regions.insert(std::make_pair(words[0], reg));
		}
		else
		{
			for (int i = stoi(words[1]); i <= stoi(words[2]); i++)
				regions[words[0]].push_back(i);
		}
	}
	ifs.close();

// read mc in the region
	int mc_crn_size = mc_crn.size(); // judge if mc_crn has established
	Pocket sca;
	sca.readrefpdb(filename);
	auto residues = sca.residues();
	std::vector<std::vector<XYZ>> ncacs;
	for (int i = 0; i < residues.size(); i++)
	{
		auto r = residues[i];
		auto reg = regions[r.chainid];
		if (auto it = find(reg.begin(), reg.end(), r.resid) == reg.end()) continue;
		std::vector<XYZ> ncac;
		ncac.push_back(getatmcrd(r, "N"));
		ncac.push_back(getatmcrd(r, "CA"));
		ncac.push_back(getatmcrd(r, "C"));
//		ncac.push_back(getatmcrd(r, "O")); // as Ocrds is arbitrary for psi
		ncacs.push_back(ncac);
		if (mc_crn_size == 0)
			mc_crn.push_back(r.chainid+"_"+to_string(r.resid)+"_"+r.resname);
	}
	mc_crds.insert(std::make_pair(filename, ncacs));
}
// read mainchain crds & lig_crds and related filenames. (for pocket)
void Theozyme::readallinfors(string filename, std::map<std::string, vector<vector<XYZ>>> &mc_crds,
		vector<string> &mc_crn, std::map<std::string, vector<vector<XYZ>>> &lig_crds, vector<string> &lig_crn)
{
// fetch all mc
	Pocket poc;
	poc.readrefpdb(filename);

	auto residues = poc.residues();
	std::vector<std::vector<XYZ>> ncacs;
	for (int i = 0; i < residues.size(); i++)
	{
		auto r = residues[i];
		std::vector<XYZ> ncac;
		ncac.push_back(getatmcrd(r, "N"));
		ncac.push_back(getatmcrd(r, "CA"));
		ncac.push_back(getatmcrd(r, "C"));
//		ncac.push_back(getatmcrd(r, "O")); // as Ocrds is arbitrary for psi
		ncacs.push_back(ncac);
		if (mc_crn.size() < residues.size())
			mc_crn.push_back(r.chainid+"_"+to_string(r.resid)+"_"+r.resname);
		else
//			assert(mc_crn[i] == r.chainid+"_"+to_string(r.resid)+"_"+r.resname);
			if (mc_crn[i] != r.chainid+"_"+to_string(r.resid)+"_"+r.resname)
			{
				std::cout << "please note that Pocket and Scaffold should only have one type with"
						" different confs. respectively." << std::endl;
				exit(1);
			}
	}
	mc_crds.insert(std::make_pair(filename, ncacs));

	auto ligands = poc.ligands();
	std::vector<std::vector<XYZ>> lcs;
	for (int i = 0; i < ligands.size(); i++)
	{
		auto l = ligands[i];
		lcs.push_back(l.atmcrds);
		if (lig_crn.size() < ligands.size())
			lig_crn.push_back(l.chainid+"_"+to_string(l.resid)+"_"+l.resname);
		else
			assert(lig_crn[i] == l.chainid+"_"+to_string(l.resid)+"_"+l.resname);
	}
	lig_crds.insert(std::make_pair(filename, lcs));
}

// see if acrds can fit with bcrds within RMSD_th
bool Theozyme::hasedge(vector<XYZ> acrds, vector<XYZ> bcrds, double RMSD_th)
{
	// acrds: pi,sj. bcrds: pa,sb.
	QuatFit qf;
	if (qf.setup(acrds, bcrds) > RMSD_th)
		return false;
	else
		return true;
}

// CMakeList has some problem that cannot define graph_new and clique_find_all.
// use clique.h to find all max_cliques in the graph.
/*
extern "C"
{
static boolean fetchsubgraph(set_t s, graph_t *g, clique_options *opts)
{
	std::vector<std::vector<int>> *members = (std::vector<std::vector<int>> *) opts->user_data;
	int i = -1;
	std::vector<int> mem;
	while ((i = set_return_next(s, i)) >= 0)
		mem.push_back(i);
	members->push_back(mem);
	return true;
}
}

vector<MatchInfo> Theozyme::graph2maxclique(vector<vector<XYZ>> pcrds, vector<vector<XYZ>> scrds, vector<double> RMSD_ths)
{
	int gsize = pcrds.size() * scrds.size();
	graph_t *g = graph_new(gsize);
	for (int i = 0; i < pcrds.size() * scrds.size() -1; i++)
		for (int j = i + 1; j < pcrds.size() * scrds.size(); j++)
		{
			int pi = i/scrds.size();
			int pj = j/scrds.size();
			int si = i%scrds.size();
			int sj = j%scrds.size();
			if (pi == pj || si == sj) continue; // one p_res can only match to one s_res
			vector<XYZ> acrds = pcrds[pi];
			acrds.insert(acrds.end(), pcrds[pj].begin(), pcrds[pj].end());
			vector<XYZ> bcrds = scrds[si];
			bcrds.insert(bcrds.end(), scrds[sj].begin(), scrds[sj].end());
			double maxr = RMSD_ths[pi] > RMSD_ths[pj] ? RMSD_ths[pi] : RMSD_ths[pj]; // only strict for two-krs
			if (hasedge(acrds, bcrds, maxr))
				GRAPH_ADD_EDGE(g, i, j);
		}
	std::vector<std::vector<int>> members;
	clique_default_options->user_function = fetchsubgraph;
	clique_default_options->user_data = (void *) (&members);
	clique_default_options->output=stderr;
	int csize = pcrds.size();
	clique_find_all(g, 1, csize, true, NULL); // max_size = pcrds.size()
	std::vector<MatchInfo> ms;
	for (auto mem : members)
	{
		if (mem.size() < pcrds.size()) continue;
		vector<int> p_s(pcrds.size(), -1); // p in s' pos
		vector<double> rmsds(pcrds.size(), 100.0);
		for (auto m : mem)
			p_s[m / scrds.size()] = m % scrds.size();
		vector<XYZ> pcrds_m;
		vector<XYZ> scrds_m;
		for (int i = 0; i < p_s.size(); i++)
		{
			if (p_s[i] != -1)
			{
				pcrds_m.insert(pcrds_m.end(), pcrds[i].begin(), pcrds[i].end());
				scrds_m.insert(scrds_m.end(), scrds[p_s[i]].begin(), scrds[p_s[i]].end());
			}
		}
		QuatFit qf;
		qf.setup(scrds_m, pcrds_m);
		MatchInfo mi;
		mi.matchposes = p_s;
		for (int i = 0; i < p_s.size(); i++)
			if (p_s[i] != -1)
			{
				auto c = pcrds[i];
				qf.transform(c);
				rmsds[i] = rmsd(c, scrds[p_s[i]]);
			}
		bool good = true;
		for (int i = 0; i < rmsds.size(); i++)
			if (rmsds[i] > RMSD_ths[i]) // filter again
				good = false;
		if (!good) continue;
		mi.qf = qf;
		mi.rmsds = rmsds;
		ms.push_back(mi);
	}

	graph_free(g);
	return ms;
}
*/
void Theozyme::graphs2mis(std::map<double, MatchInfo> &mis, vector<MatchInfo> &ms,
		std::string pfname, std::string sfname, vector<double> RMSD_ths, set<int> kr_ids)
{
	for (auto m : ms)
	{
		double ave_r = 0.0;
		for (int i = 0; i < m.rmsds.size(); i++)
		{
			double w = 1;
			if (kr_ids.count(i) > 0) w = 100; // assume KR_ratio = 0.1, OR_ratio = 1, then *100 makes 1KR ~ 10OR.
			double ratio = m.rmsds[i]/RMSD_ths[i];
			if (ratio > 1)
				ave_r += w*5;
			else
				ave_r += w*(1-cos(M_PI*ratio))/2;
		}
		m.pocket_filename = pfname;
		m.scaffold_filename = sfname;
		mis.insert(std::make_pair(ave_r/RMSD_ths.size(), m)); // Wkr:Wor = 1000:1.
	}
}
/*
// pocket can randomly move during match
void Theozyme::graphicmatch_ka(vector<KeyAtom> poc_ka, vector<KeyAtom> sca_ka, vector<double> RMSD_ths,
		std::map<double, MatchInfo> &mis, string pfilename, string sfilename)
{
	std::map<int, pair<int, int>> v2p; // vid -> poc, sca. Used for minimize the graph.
	int k = 0;
	for (int i = 0; i < poc_ka.size(); i++)
		for (int j = 0; j < sca_ka.size(); j++)
			if (poc_ka[i].rname == sca_ka[j].rname && poc_ka[i].id == sca_ka[j].id)
			{
				pair<int, int> ps;
				ps.first = i;
				ps.second = j;
				v2p.insert(make_pair(k, ps));
				k++;
			}
	int gsize = k;
	graph_t *g = graph_new(gsize);
	for (int i = 0; i < k-1; i++)
		for (int j = i+1; j < k; j++)
		{
			int pi = v2p[i].first;
			int si = v2p[i].second;
			int pj = v2p[j].first;
			int sj = v2p[j].second;
			// one p_res can only match to (non-overlapping) one s_res
			if (pi == pj) continue;
			if (sca_ka[si].cid == sca_ka[sj].cid && sca_ka[si].rid == sca_ka[sj].rid) continue;
			vector<XYZ> acrds = poc_ka[pi].crds;
			acrds.insert(acrds.end(), poc_ka[pj].crds.begin(), poc_ka[pj].crds.end());
			vector<XYZ> bcrds = sca_ka[si].crds;
			bcrds.insert(bcrds.end(), sca_ka[sj].crds.begin(), sca_ka[sj].crds.end());
			double maxr = RMSD_ths[pi] > RMSD_ths[pj] ? RMSD_ths[pi] : RMSD_ths[pj]; // only strict for two-krs
			if (hasedge(acrds, bcrds, maxr))
				GRAPH_ADD_EDGE(g, i, j);
		}
	std::vector<std::vector<int>> members;
	clique_default_options->user_function = fetchsubgraph;
	clique_default_options->user_data = (void *) (&members);
	clique_default_options->output=stderr;
	int csize = poc_ka.size();
	clique_find_all(g, csize, csize, true, NULL);
	// max_size = pcrds.size(), 1->csize s.t. finding matches for all pocket residues.
	std::vector<MatchInfo> ms;
	for (auto mem : members)
	{
		if (mem.size() < poc_ka.size())
			continue;
		vector<int> p_s(poc_ka.size(), -1); // p in s' pos
		vector<double> rmsds(poc_ka.size(), 100.0);
		for (auto m : mem)
			p_s[v2p[m].first] = v2p[m].second;
		vector<XYZ> pcrds_m;
		vector<XYZ> scrds_m;
		for (int i = 0; i < p_s.size(); i++)
		{
			if (p_s[i] != -1)
			{
				pcrds_m.insert(pcrds_m.end(), poc_ka[i].crds.begin(), poc_ka[i].crds.end());
				scrds_m.insert(scrds_m.end(), sca_ka[p_s[i]].crds.begin(), sca_ka[p_s[i]].crds.end());
			}
		}
		QuatFit qf;
		qf.setup(scrds_m, pcrds_m);
		MatchInfo mi;
		mi.matchposes = p_s;
		for (int i = 0; i < p_s.size(); i++)
			if (p_s[i] != -1)
			{
				vector<XYZ> c = poc_ka[i].crds;
				qf.transform(c);
				rmsds[i] = rmsd(c, sca_ka[p_s[i]].crds);
			}
		bool good = true;
		for (int i = 0; i < rmsds.size(); i++)
			if (rmsds[i] > RMSD_ths[i]) // filter again
				good = false;
		if (!good) continue;
		mi.qf = qf;
		mi.rmsds = rmsds;
		ms.push_back(mi);
	}
	graph_free(g);
	set<int> kr_ids;
	graphs2mis(mis, ms, pfilename, sfilename, RMSD_ths, kr_ids);
}
*/
// kr using graphicmatch, while others using directmatch.
void Theozyme::graphicmatch_ka(vector<KeyAtom> poc_ka, vector<KeyAtom> sca_ka,
		vector<double> RMSD_ths, std::map<double, MatchInfo> &mis, string pfilename, string sfilename,
		set<int> kr_ids, int keyresredconf)
{
	std::map<int, pair<int, int>> v2p; // vid -> poc, sca. Used for minimize the graph.
	int k = 0;
	for (int i = 0; i < poc_ka.size(); i++)
	{
		if (kr_ids.count(i) == 0) continue;
		for (int j = 0; j < sca_ka.size(); j++)
			if (poc_ka[i].rname == sca_ka[j].rname && poc_ka[i].id == sca_ka[j].id)
			{
				pair<int, int> ps;
				ps.first = i;
				ps.second = j;
				v2p.insert(make_pair(k, ps));
				k++;
			}
	}
	int gsize = k;
	Graph g;
	g.graph_new(gsize);
	for (int i = 0; i < k-1; i++)
		for (int j = i+1; j < k; j++)
		{
			int pi = v2p[i].first;
			int si = v2p[i].second;
			int pj = v2p[j].first;
			int sj = v2p[j].second;
			// one p_res can only match to (non-overlapping) one s_res
			if (pi == pj) continue;
			if (sca_ka[si].cid == sca_ka[sj].cid && sca_ka[si].rid == sca_ka[sj].rid) continue;
			vector<XYZ> acrds = poc_ka[pi].crds;
			acrds.insert(acrds.end(), poc_ka[pj].crds.begin(), poc_ka[pj].crds.end());
			vector<XYZ> bcrds = sca_ka[si].crds;
			bcrds.insert(bcrds.end(), sca_ka[sj].crds.begin(), sca_ka[sj].crds.end());
			double maxr = RMSD_ths[pi] > RMSD_ths[pj] ? RMSD_ths[pi] : RMSD_ths[pj]; // only strict for two-krs
			if (hasedge(acrds, bcrds, maxr))
				g.add_edge(i, j);
		}
	std::vector<std::vector<int>> members;
	vector<MyClique::Clique> cs = findclique_map(g, kr_ids.size());
	g.graph_free();
	for (auto c : cs)
	{
		vector<int> member;
		for (auto i : c.ids)
			member.push_back(i);
		members.push_back(member);
	}
	cs.clear();
	// judge if non-redundant
	std::map<double, vector<int>> krmsd_memid; // kr's rmsd ranking for each member.
	for (int n = 0; n < members.size(); n++)
	{
		auto mem = members[n];
		vector<XYZ> pcrds_m, scrds_m;
		for (auto m : mem)
		{
			pcrds_m.insert(pcrds_m.end(), poc_ka[v2p[m].first].crds.begin(), poc_ka[v2p[m].first].crds.end());
			scrds_m.insert(scrds_m.end(), sca_ka[v2p[m].second].crds.begin(), sca_ka[v2p[m].second].crds.end());
		}
		QuatFit qf;
		qf.setup(scrds_m, pcrds_m);

		// requiring kr must within RMSD_ths
		bool within = true;
		for (auto m : mem)
		{
			vector<XYZ> kr_crd = poc_ka[v2p[m].first].crds;
			qf.transform(kr_crd);
			if (rmsd(kr_crd, sca_ka[v2p[m].second].crds) >= RMSD_ths[v2p[m].first])
			{
				within = false;
				break;
			}
		}
		if (!within) continue;

		qf.transform(pcrds_m);
		double r = rmsd(scrds_m, pcrds_m);
		if (krmsd_memid.count(r) > 0)
			krmsd_memid[r].push_back(n);
		else
		{
			vector<int> new_n {n};
			krmsd_memid.insert(make_pair(r, new_n));
		}
	}
	cout << "Initially there are " << members.size() << " key residue matches." << endl;
	vector<vector<string>> nr_codes;
	vector<int> nr_codes_conf; // for each code, confs could be no more than keyresredconf.
	vector<vector<int>> members_nr; // non-redundant members
	for (auto it = krmsd_memid.begin(); it != krmsd_memid.end(); it++)
	{
		for (auto mid : it->second)
		{
			bool add_c = true;
			vector<string> new_code;
			for (auto m : members[mid])
			{
				int sid = v2p[m].second;
				int pid = v2p[m].first;
				new_code.push_back(to_string(pid) + ":" + sca_ka[sid].cid + "_" + to_string(sca_ka[sid].rid));
			}
			if (nr_codes.size() == 0)
			{
				nr_codes.push_back(new_code);
				nr_codes_conf.push_back(1);
			}
			else
			{
				int n_c = 0;
				for (auto nrc : nr_codes)
				{
					int same = 0;
					for (int n = 0; n < new_code.size(); n++)
						if (nrc[n] == new_code[n]) same++;
					if (same == new_code.size()) n_c++;
					if (n_c >= keyresredconf)
					{
						add_c = false;
						break;
					}
				}
			}
			if (add_c)
			{
				nr_codes.push_back(new_code);
				members_nr.push_back(members[mid]);
			}
		}
	}
	cout << "There are " << members_nr.size() << " non-redundant matches for key residues." << endl;
	std::vector<MatchInfo> ms;
	vector<vector<string>> nr_members; // non_redundant members.
	for (auto mem : members_nr)
	{
		assert(mem.size() == kr_ids.size());
		vector<int> p_s(poc_ka.size(), -1); // p in s' pos
		for (auto m : mem)
			p_s[v2p[m].first] = v2p[m].second;
		vector<XYZ> pcrds_m;
		vector<XYZ> scrds_m;
		for (int i = 0; i < p_s.size(); i++)
		{
			if (p_s[i] != -1)
			{
				pcrds_m.insert(pcrds_m.end(), poc_ka[i].crds.begin(), poc_ka[i].crds.end());
				scrds_m.insert(scrds_m.end(), sca_ka[p_s[i]].crds.begin(), sca_ka[p_s[i]].crds.end());
			}
		}
		QuatFit qf; // directly using key resiudes' qf for the whole pocket.
		qf.setup(scrds_m, pcrds_m);
		// update new poc_kas for this match(qf).
		vector<KeyAtom> poc_ka_new; // = poc_ka * qf.
		for (auto pk : poc_ka)
			poc_ka_new.push_back(pk);
		for (auto& pkn : poc_ka_new)
			qf.transform(pkn.crds);
		// find positions for other residues.
		vector<vector<int>> p_ss = directmatch_oka(poc_ka_new, sca_ka, RMSD_ths, p_s);
		for (auto pss : p_ss)
		{
			vector<double> rmsds(RMSD_ths.size(), 100);
			MatchInfo mi;
			mi.matchposes = pss;
			for (int i = 0; i < pss.size(); i++)
				if (pss[i] != -1)
				{
					vector<XYZ> c = poc_ka[i].crds;
					qf.transform(c);
					rmsds[i] = rmsd(c, sca_ka[pss[i]].crds);
				}
			mi.qf = qf;
			mi.rmsds = rmsds;
			ms.push_back(mi);
		}
	}

	graphs2mis(mis, ms, pfilename, sfilename, RMSD_ths, kr_ids);
}

/*
// kr using graphicmatch-like method without qf, all using directmatch.
void Theozyme::directmatch_ka_graph(vector<KeyAtom> poc_ka, vector<KeyAtom> sca_ka,
		vector<double> RMSD_ths, std::map<double, MatchInfo> &mis, string pfilename, string sfilename,
		set<int> kr_ids)
{
	std::map<int, pair<int, int>> v2p; // vid -> poc, sca. Used for minimize the graph.
	int k = 0;
	for (int i = 0; i < poc_ka.size(); i++)
	{
		if (kr_ids.count(i) == 0) continue;
		for (int j = 0; j < sca_ka.size(); j++)
			if (poc_ka[i].rname == sca_ka[j].rname && poc_ka[i].id == sca_ka[j].id)
			{
				pair<int, int> ps;
				ps.first = i;
				ps.second = j;
				v2p.insert(make_pair(k, ps));
				k++;
			}
	}
	int gsize = k;
	graph_t *g = graph_new(gsize);
	for (int i = 0; i < k-1; i++)
		for (int j = i+1; j < k; j++)
		{
			int pi = v2p[i].first;
			int si = v2p[i].second;
			int pj = v2p[j].first;
			int sj = v2p[j].second;
			// one p_res can only match to (non-overlapping) one s_res
			if (pi == pj) continue;
			if (sca_ka[si].cid == sca_ka[sj].cid && sca_ka[si].rid == sca_ka[sj].rid) continue;
			vector<XYZ> acrds = poc_ka[pi].crds;
			acrds.insert(acrds.end(), poc_ka[pj].crds.begin(), poc_ka[pj].crds.end());
			vector<XYZ> bcrds = sca_ka[si].crds;
			bcrds.insert(bcrds.end(), sca_ka[sj].crds.begin(), sca_ka[sj].crds.end());
			double maxr = RMSD_ths[pi] > RMSD_ths[pj] ? RMSD_ths[pi] : RMSD_ths[pj]; // only strict for two-krs
			if (rmsd(acrds, bcrds) > maxr)
				GRAPH_ADD_EDGE(g, i, j);
		}
	std::vector<std::vector<int>> members;
	clique_default_options->user_function = fetchsubgraph;
	clique_default_options->user_data = (void *) (&members);
	clique_default_options->output=stderr;
	int csize = kr_ids.size();
	clique_find_all(g, csize, csize, true, NULL);
	graph_free(g);

	// judge if non-redundant
	std::map<double, vector<int>> krmsd_memid; // kr's rmsd ranking for each member.
	for (int n = 0; n < members.size(); n++)
	{
		auto mem = members[n];
		vector<XYZ> pcrds_m, scrds_m;
		for (auto m : mem)
		{
			pcrds_m.insert(pcrds_m.end(), poc_ka[v2p[m].first].crds.begin(), poc_ka[v2p[m].first].crds.end());
			scrds_m.insert(scrds_m.end(), sca_ka[v2p[m].second].crds.begin(), sca_ka[v2p[m].second].crds.end());
		}
		double r = rmsd(scrds_m, pcrds_m);
		if (krmsd_memid.count(r) > 0)
			krmsd_memid[r].push_back(n);
		else
		{
			vector<int> new_n {n};
			krmsd_memid.insert(make_pair(r, new_n));
		}
	}
	cout << "Initially there are " << members.size() << " key residue matches." << endl;
	vector<vector<string>> nr_codes;
	vector<vector<int>> members_nr; // non-redundant members
	for (auto it = krmsd_memid.begin(); it != krmsd_memid.end(); it++)
	{
		for (auto mid : it->second)
		{
			bool new_c = true;
			vector<string> new_code;
			for (auto m : members[mid])
			{
				int sid = v2p[m].second;
				new_code.push_back(sca_ka[sid].cid + "_" + to_string(sca_ka[sid].rid));
			}
			if (nr_codes.size() == 0)
				nr_codes.push_back(new_code);
			else
				for (auto nrc : nr_codes)
				{
					int same = 0;
					for (int n = 0; n < new_code.size(); n++)
						if (nrc[n] == new_code[n]) same++;
					if (same == new_code.size()) new_c = false;
					if (!new_c) break;
				}
			if (new_c)
			{
				nr_codes.push_back(new_code);
				members_nr.push_back(members[mid]);
			}
		}
	}
	cout << "There are " << members_nr.size() << " non-redundant matches for key residues." << endl;
	std::vector<MatchInfo> ms;
	vector<vector<string>> nr_members; // non_redundant members.
	for (auto mem : members_nr)
	{
		assert(mem.size() == kr_ids.size());
		vector<int> p_s(poc_ka.size(), -1); // p in s' pos
		for (auto m : mem)
			p_s[v2p[m].first] = v2p[m].second;
		vector<XYZ> pcrds_m;
		vector<XYZ> scrds_m;
		for (int i = 0; i < p_s.size(); i++)
		{
			if (p_s[i] != -1)
			{
				pcrds_m.insert(pcrds_m.end(), poc_ka[i].crds.begin(), poc_ka[i].crds.end());
				scrds_m.insert(scrds_m.end(), sca_ka[p_s[i]].crds.begin(), sca_ka[p_s[i]].crds.end());
			}
		}
		QuatFit qf; // no move.
		// find positions for other residues.
		vector<vector<int>> p_ss = directmatch_oka(poc_ka, sca_ka, RMSD_ths, p_s);
		for (auto pss : p_ss)
		{
			vector<double> rmsds(poc_ka.size(), 100.0);
			MatchInfo mi;
			mi.matchposes = pss;
			for (int i = 0; i < pss.size(); i++)
				if (pss[i] != -1)
				{
					vector<XYZ> c = poc_ka[i].crds;
					rmsds[i] = rmsd(c, sca_ka[pss[i]].crds);
				}
			mi.qf = qf;
			mi.rmsds = rmsds;
			ms.push_back(mi);
		}
	}
	graphs2mis(mis, ms, pfilename, sfilename, RMSD_ths, kr_ids);
}
*/
void Theozyme::addmatchinfo(vector<MatchInfo> &ms, vector<KeyAtom> sca_ka, vector<vector<int>> pidsid, vector<vector<double>> ps_rmsd,
		int pid, set<int> kr_ids, vector<double> RMSD_ths)
{
	int kr_num = pidsid[pid].size();
	if (ms.size() == 0)
	{
		while(kr_num > 0)
		{
			MatchInfo m;
			int ith = pidsid[pid].size() - kr_num;
			string scode = sca_ka[pidsid[pid][ith]].cid + "_" + to_string(sca_ka[pidsid[pid][ith]].rid);
			m.scacodes.insert(scode);
			vector<int> matchposes(pidsid.size(), -1);
			vector<double> rmsds(RMSD_ths.size(), 100);
			matchposes[pid] = pidsid[pid][ith];
			rmsds[pid] = ps_rmsd[pid][ith];
			m.matchposes.assign(matchposes.begin(), matchposes.end());
			m.rmsds.assign(rmsds.begin(), rmsds.end());
			ms.push_back(m);
			kr_num--;
		}
		if (kr_ids.count(pid) == 0)
		{
			MatchInfo m;
			vector<int> matchposes(pidsid.size(), -1);
			vector<double> rmsds(RMSD_ths.size(), 100);
			m.matchposes.assign(matchposes.begin(), matchposes.end());
			m.rmsds.assign(rmsds.begin(), rmsds.end());
			ms.push_back(m);
		}
	}
	else
	{
		vector<MatchInfo> ms_add;
		for (auto m : ms)
		{
			if (kr_ids.count(pid) == 0)
			{
				ms_add.push_back(m);
			}
			while(kr_num > 0)
			{
				int ith = pidsid[pid].size() - kr_num;
				string scode = sca_ka[pidsid[pid][ith]].cid + "_" + to_string(sca_ka[pidsid[pid][ith]].rid);
				if (m.scacodes.count(scode) > 0)
				{
					kr_num--;
					continue;
				}
				else
				{
					MatchInfo m_add;
					m_add.matchposes.assign(m.matchposes.begin(), m.matchposes.end());
					m_add.rmsds.assign(m.rmsds.begin(), m.rmsds.end());
					for (auto sc : m.scacodes)
						m_add.scacodes.insert(sc);
					m_add.matchposes[pid] = pidsid[pid][ith];
					m_add.rmsds[pid] = ps_rmsd[pid][ith];
					m_add.scacodes.insert(scode);
					ms_add.push_back(m_add);
					kr_num--;
				}
			}
			kr_num = pidsid[pid].size(); // for the next m, kr_num should go over again.
		} // each old m.
		ms.clear();
		for (auto m_add : ms_add)
			ms.push_back(m_add);
	}
}


// all using directmatch
void Theozyme::directmatch_ka(vector<KeyAtom> poc_ka, vector<KeyAtom> sca_ka,
		vector<double> RMSD_ths, std::map<double, MatchInfo> &mis, string pfilename, string sfilename,
		set<int> kr_ids, int keyresredconf)
{
	vector<vector<int>> pidsid(poc_ka.size()); // <every_pid<possible_sid>>
	vector<vector<double>> ps_rmsd(poc_ka.size()); // <every_pid<rmsd>>
	// in each Scaffild_site, the key residue should only keep the MinRMSD match result for acceleration.
	bool kr_find = true; // all key residues must find match.

	for (auto ki : kr_ids)
	{
		map<string, map<double, int>> site_krconf; // for each site, kr's conf number.
		double minRMSD = 100.0;
		auto pk = poc_ka[ki];
		for (int j = 0; j < sca_ka.size(); j++)
		{
			auto sk = sca_ka[j];
			if (pk.rname == sk.rname && pk.id == sk.id)
			{
				double r = rmsd(pk.crds, sk.crds);
				if (r < RMSD_ths[ki])
				{
					string site = sk.cid + "_" + to_string(sk.rid);
					if (site_krconf.count(site) > 0)
					{
						auto &krconf = site_krconf[site];
						if (krconf.size() < keyresredconf)
							krconf.insert(make_pair(r, j));
						else
						{
							auto iter = krconf.end();
							iter--;
							if (iter->first > r)
							{
								krconf.erase(iter);
								krconf.insert(make_pair(r, j));
							} // always keep the top keyresredconf. confs.
						}
					}
					else
					{
						map<double, int> krconf;
						krconf.insert(make_pair(r, j));
						site_krconf.insert(make_pair(site, krconf));
					}
				}
			}
		}
		if (site_krconf.size() == 0)
		{
			cout << "Key Residue " << ki << " didn't find match (counted from 0)." << endl;
			kr_find = false;
		}
		else
		{
			cout << "Key Residue " << ki << " find sites: ";
			for (auto ita = site_krconf.begin(); ita != site_krconf.end(); ita++)
			{
				auto krconf = ita->second;
				for (auto iter = krconf.begin(); iter != krconf.end(); iter++)
				{
					pidsid[ki].push_back(iter->second);
					ps_rmsd[ki].push_back(iter->first);
					cout << sca_ka[iter->second].cid << " " << sca_ka[iter->second].rid << "; ";
				}
			}
			cout << endl;
		}
	}
	if (!kr_find) exit(1);
	for (int i = 0; i < poc_ka.size(); i++)
	{
		if (kr_ids.count(i) > 0) continue;
		for (int j = 0; j < sca_ka.size(); j++)
		{
			auto pk = poc_ka[i];
			auto sk = sca_ka[j];
			if (pk.rname == sk.rname && pk.id == sk.id)
			{
				double r = rmsd(pk.crds, sk.crds);
				if (r < RMSD_ths[i])
				{
					pidsid[i].push_back(j);
					ps_rmsd[i].push_back(r);
				}
			}
		}
	}

	// new version.
	vector<MatchInfo> ms; // final matchinfo for mis.
	for (auto ki : kr_ids)
		addmatchinfo(ms, sca_ka, pidsid, ps_rmsd, ki, kr_ids, RMSD_ths);
	for (int i = 0; i < pidsid.size(); i++)
		if (kr_ids.count(i) == 0)
			addmatchinfo(ms, sca_ka, pidsid, ps_rmsd, i, kr_ids, RMSD_ths);
/*
 // old version.
	vector<MatchInfo> ms; // final matchinfo for mis.
	if (pidsid[0].size() == 0)
	{
		MatchInfo m;
		m.matchposes.push_back(-1);
		m.rmsds.push_back(10.0);
		ms.push_back(m);
		cout << " Pocket Residue 0 doesn't find matche." << endl;
	}
	else
	{
		for (int si = 0; si < pidsid[0].size(); si++)
		{
			MatchInfo m;
			m.matchposes.push_back(pidsid[0][si]);
			m.rmsds.push_back(ps_rmsd[0][si]);
			ms.push_back(m);
		}
		cout << " Pocket Residue 0 has " << pidsid[0].size() << " matches." << endl;
	}
	vector<MatchInfo> ms_temp; // temperory ms.
	int pi = 1;
	while(true)
	{
		if (pidsid[pi].size() == 0)
		{
			cout << " Pocket Residue " << pi << " doesn't find matche." << endl;
			for (auto &m : ms)
			{
				m.matchposes.push_back(-1);
				m.rmsds.push_back(10.0);
			}
		}
		else
		{
			cout << " Pocket Residue " << pi << " has " << pidsid[pi].size() << " matches." << endl;
			int msize = ms.size();
			if (kr_ids.count(pi) > 0)
			{
				for (int i = 1; i < pidsid[pi].size(); i++)
					for (int j = 0; j < msize; j++)
						ms.push_back(ms[j]);
				for (int i = 0; i < pidsid[pi].size(); i++)
					for (int j = 0; j < msize; j++)
					{
						ms[i*msize+j].matchposes.push_back(pidsid[pi][i]);
						ms[i*msize+j].rmsds.push_back(ps_rmsd[pi][i]);
					}
			}
			else // other residue can have no match.
			{
				for (int i = 0; i < pidsid[pi].size(); i++)
					// first is for -1 (didn't find match), so it need one more copy.
					for (int j = 0; j < msize; j++)
						ms.push_back(ms[j]);
				for (int i = 0; i < pidsid[pi].size(); i++)
					for (int j = 0; j < msize; j++)
					{
						ms[i*msize+j].matchposes.push_back(pidsid[pi][i]);
						ms[i*msize+j].rmsds.push_back(ps_rmsd[pi][i]);
					}
				for (int j = 0; j < msize; j++)
				{
					ms[pidsid[pi].size()*msize+j].matchposes.push_back(-1);
					ms[pidsid[pi].size()*msize+j].rmsds.push_back(10.0);
				} // for -1.
			}
		}
		pi++;
		if (pi == pidsid.size()) break;
	}
*/
	for (auto &mss : ms)
		mss.dm = true;

	graphs2mis(mis, ms, pfilename, sfilename, RMSD_ths, kr_ids);
}

vector<vector<int>> Theozyme::directmatch_oka(vector<KeyAtom> poc_ka_new, vector<KeyAtom> sca_ka, vector<double> RMSD_ths,
		vector<int> p_s)
{
	std::map<int, vector<int>> p_ss; // recording possible sid in each pid.
	set<string> existings;
	for (auto s : p_s)
		if (s != -1)
			existings.insert(sca_ka[s].cid + "_" + to_string(sca_ka[s].rid));
	for (int p = 0; p < p_s.size(); p++)
	{
		int s = p_s[p];
		if (s != -1)
		{
			vector<int> ss {s};
			p_ss.insert(make_pair(p, ss));
		}
		else
		{
			vector<int> ss {-1};
			auto pkn = poc_ka_new[p];
			for (int a = 0; a < sca_ka.size(); a++)
			{
				auto sk = sca_ka[a];
				string skcode = sk.cid+"_"+to_string(sk.rid);
				if (existings.count(skcode) > 0) continue;
				if (pkn.rname == sk.rname && pkn.id == sk.id)
					if (rmsd(pkn.crds, sk.crds) < RMSD_ths[p])
						ss.push_back(a);
			}
			p_ss.insert(make_pair(p, ss));
		}
	}
//
//	// test
//	for (int i = 0; i < p_ss.size(); i++)
//	{
//		cout << i << " : ";
//		for (auto p : p_ss[i])
//			if (p == -1)
//				cout << "Nan; ";
//			else
//				cout << sca_ka[p].rid << "; ";
//		cout << endl;
//	}

	vector<int> first_p_s_g(p_s.size(), -1);
	for (int p = 0; p < p_s.size(); p++)
		if (p_s[p] != -1)
			first_p_s_g[p] = p_s[p]; // first do KR
	vector<vector<int>> p_s_grps{first_p_s_g}; // one p_s for kr can derive a p_s_grps.
	for (int p = 0; p < p_s.size(); p++)
		if (p_s[p] == -1)
			addmi_simple(p_s_grps, p, p_ss[p], sca_ka, p_s.size()); // then do OR

	return p_s_grps;
}

void Theozyme::addmi_simple(vector<vector<int>> &p_s_grps, int p, vector<int> sites, vector<KeyAtom> sca_ka, int p_num)
{
	set<string> exists;
	vector<vector<int>> p_s_grps_new;
	if (p_s_grps.size() == 0)
		for (auto s : sites)
		{
			vector<int> psg(p_num, -1);
			psg[p] == s;
			p_s_grps_new.push_back(psg);
		}
	else
		for (auto psg : p_s_grps)
			for (auto s : sites)
			{
				set<string> exists;
				for (auto i : psg)
					if (i != -1)
					{
						string si = sca_ka[i].cid + "_" + to_string(sca_ka[i].rid);
						exists.insert(si);
					}
				if (s != -1 && exists.count(sca_ka[s].cid+"_"+to_string(sca_ka[s].rid)) > 0) continue;
				psg[p] = s;
				p_s_grps_new.push_back(psg);
			}
	p_s_grps.clear();
	for (auto grp : p_s_grps_new)
		p_s_grps.push_back(grp);
}

// branch and bound method for ka match
void Theozyme::bnb_ka(vector<KeyAtom> poc_ka, vector<KeyAtom> sca_ka,
		vector<double> RMSD_ths, std::map<double, MatchInfo> &mis, string pfilename, string sfilename)
{
	// sort from small to large
	std::map<int, vector<int>> cn_pids; // candi number (sort) with its poc_id
	vector<vector<int>> candis; // pid -> candi sid
	vector<int> candi_nums(poc_ka.size(), 0); // candi numbers
	for (int i = 0; i < poc_ka.size(); i++)
	{
		vector<int> candi;
		for (int j = 0; j < sca_ka.size(); j++)
			if (poc_ka[i].rname == sca_ka[j].rname && poc_ka[i].id == sca_ka[j].id)
			{
				candi.push_back(j);
				candi_nums[i]++;
			}
		candis.push_back(candi);
		if (cn_pids.count(candi_nums[i]) > 0)
			cn_pids[candi_nums[i]].push_back(i);
		else
		{
			vector<int> pids{i};
			cn_pids.insert(make_pair(candi_nums[i], pids));
		}
	}
	vector<int> sort_pid; // cn_pids->sort_pid for convenience.
	for (auto it = cn_pids.begin(); it != cn_pids.end(); it++)
		for (auto c : it->second)
			sort_pid.push_back(c);

	// establish edge map for convenient judgement
	std::map<string, set<string>> ps_map; // {pi,sa}<-edge->{pj,sb}
	for (int i = 0; i < candis.size()-1; i++)
		for (int j = i+1; j < candis.size(); j++)
			for (auto a : candis[i])
				for (auto b : candis[j])
				{
					if (sca_ka[a].cid == sca_ka[b].cid && sca_ka[a].rid == sca_ka[b].rid) continue;
					vector<XYZ> acrds = poc_ka[i].crds;
					acrds.insert(acrds.end(), poc_ka[j].crds.begin(), poc_ka[j].crds.end());
					vector<XYZ> bcrds = sca_ka[a].crds;
					bcrds.insert(bcrds.end(), sca_ka[b].crds.begin(), sca_ka[b].crds.end());
					double maxr = RMSD_ths[i] > RMSD_ths[j] ? RMSD_ths[i] : RMSD_ths[j]; // only strict for two-krs
					if (hasedge(acrds, bcrds, maxr))
					{
						string ia = to_string(i) + "_" + to_string(a);
						string jb = to_string(j) + "_" + to_string(b);
						if (ps_map.count(ia) > 0)
							ps_map[ia].insert(jb);
						else
						{
							set<string> ias;
							ias.insert(jb);
							ps_map.insert(make_pair(ia, ias));
						}
						if (ps_map.count(jb) > 0)
							ps_map[jb].insert(ia);
						else
						{
							set<string> jbs;
							jbs.insert(ia);
							ps_map.insert(make_pair(jb, jbs));
						}
					}
				}

	// branch and bound method
	int po = 0; // dealing with po and its sons pi.
	vector<vector<int>> member;
	for (auto ca : candis[sort_pid[po]])
	{
		vector<int> mem {ca};
		member.push_back(mem);
	} // initialize by first.
	cout << "Initial " << po << " with size " << member.size() << endl;
	while (po < poc_ka.size()-1 && member.size() > 0) // Try adding po+1 into previous maximal clique
	{
		vector<vector<int>> mem_new;
		int pi = po + 1;
		cout << "Try adding " << pi << " with size " << candis[sort_pid[pi]].size()
				<< " into " << member.size() << endl;
		for (auto mem : member)
		{
			// expand or delete;
			for (auto sb : candis[sort_pid[pi]])
			{
				bool mc = true; // maximal_clique
				string ib = to_string(sort_pid[pi]) + "_" + to_string(sb);
				for (int o = 0; o < mem.size(); o++)
				{
					string oa = to_string(sort_pid[o]) + "_" + to_string(mem[o]);
					assert(ps_map.count(ib) >0 && ps_map.count(oa) > 0);
					if (ps_map[ib].count(oa) == 0)
					{
						mc = false;
						break;
					}
				}
				if (mc)
				{
					vector<int> mem_n;
					for (auto m : mem) mem_n.push_back(m);
					mem_n.push_back(sb);
					mem_new.push_back(mem_n);
				} // update mem_new
			}
		}
		member.clear();
		for (auto mem_n : mem_new)
			member.push_back(mem_n);
		po++;
	}

	// update mis
	std::vector<MatchInfo> ms;
	for (auto mem : member)
	{
		assert(mem.size() == poc_ka.size());
		vector<int> p_s(poc_ka.size(), -1); // p in s' pos
		vector<double> rmsds(poc_ka.size(), 100.0);
		for (int i = 0; i < mem.size(); i++)
			p_s[sort_pid[i]] = mem[i];
		vector<XYZ> pcrds_m;
		vector<XYZ> scrds_m;
		for (int i = 0; i < p_s.size(); i++)
		{
			if (p_s[i] != -1)
			{

				pcrds_m.insert(pcrds_m.end(), poc_ka[i].crds.begin(), poc_ka[i].crds.end());
				scrds_m.insert(scrds_m.end(), sca_ka[p_s[i]].crds.begin(), sca_ka[p_s[i]].crds.end());
			}
		}
		QuatFit qf;
		qf.setup(scrds_m, pcrds_m);
		MatchInfo mi;
		mi.matchposes = p_s;
		for (int i = 0; i < p_s.size(); i++)
			if (p_s[i] != -1)
			{
				vector<XYZ> c = poc_ka[i].crds;
				qf.transform(c);
				rmsds[i] = rmsd(c, sca_ka[p_s[i]].crds);
			}
		bool good = true;
		for (int i = 0; i < rmsds.size(); i++)
			if (rmsds[i] > RMSD_ths[i]) // filter again
				good = false;
		if (!good) continue;
		mi.qf = qf;
		mi.rmsds = rmsds;
		ms.push_back(mi);
	}
	set<int> kr_ids;
	graphs2mis(mis, ms, pfilename, sfilename, RMSD_ths, kr_ids);
}

std::map<double, MatchInfo> Theozyme::checkclash_lig_sc(std::map<double, MatchInfo> mis,
		std::map<string, vector<vector<XYZ>>> sca_crds_all)
{
	std::map<double, MatchInfo> mis_new;
	for (auto itm = mis.begin(); itm != mis.end(); itm++)
	{
		MatchInfo mi = itm->second;
		// sca_mc;
		string sf = mi.scaffold_filename;
		vector<vector<XYZ>> sca_mc = sca_crds_all[sf];

		// poc_lig;
		string pf = mi.pocket_filename;
		Pocket pocket;
		pocket.readrefpdb(pf);
		auto qf = itm->second.qf;
		pocket.transform(qf);
		auto poc_ligs = pocket.ligands();
		bool lig_clash = false;
		for (auto lig : poc_ligs)
		{
			for (auto la : lig.atmcrds)
			{
				for (auto sm : sca_mc)
				{
					for (auto s : sm)
					{
						if (la.distance(s) < 3.0)
							lig_clash = true;
						if (lig_clash) break;
					}
					if (lig_clash) break;
				}
				if (lig_clash) break;
			}
			if (lig_clash) break;
		}
		if (!lig_clash)	// copy mis_member
			mis_new.insert(std::make_pair(itm->first, mi)); // <score, mi>
	}
	return mis_new;
}

std::map<double, MatchInfo> Theozyme::check_write_poc(std::map<double, MatchInfo> mis_n1,
		std::map<string, vector<KeyAtom>> sca_kas)
{
	std::map<double, MatchInfo> mi_new;
	int i = 0;
	for (auto itm = mis_n1.begin(); itm != mis_n1.end(); itm++)
	{
		MatchInfo mi = itm->second;
		// do transform & write pdb.
		Pocket pocket;
		pocket.readrefpdb(mi.pocket_filename);
		auto qf = mi.qf;
		pocket.transform(qf);
		// change_sc;
		renew_residues(pocket, sca_kas[mi.scaffold_filename], mi);
		// judge clash between res-res (3 A) and lig-res (1.5 A).
		bool clash = false;
		auto residues = pocket.residues();
		auto ligands = pocket.ligands();
		auto mp = mi.matchposes;
		for (int p = 0; p < residues.size()-1; p++)
		{
			if (mp[p] < 0) continue;
			auto rp = residues[p];
			for (int q = p+1; q < residues.size(); q++)
			{
				if (mp[q] < 0) continue;
				auto rq = residues[q];
				for (auto rp_c : rp.atmcrds)
				{
					for (auto rq_c : rq.atmcrds)
						if (rp_c.distance(rq_c) < 3.0)
						{
							clash = true;
							break;
						}
					if (clash) break;
				}
				if (clash) break;
			}
			if (clash) break;
		}
		for (int p = 0; p < residues.size(); p++)
		{
			if (mp[p] < 0) continue;
			auto rp = residues[p];
			for (int k = 0; k < ligands.size(); k++)
			{
				auto lk = ligands[k];
				for (auto rp_c : rp.atmcrds)
				{
					for (auto lk_c : lk.atmcrds)
						if (rp_c.distance(lk_c) < 1.5)
						{
							clash = true;
							break;
						}
					if (clash) break;
				}
				if (clash) break;
			}
			if (clash) break;
		}
		if (!clash)
		{
			for (int k = 0; k < ligands.size(); k++)
			{
				auto lk = ligands[k];
				for (auto rp_c : residues[residues.size()-1].atmcrds)
				{
					for (auto lk_c : lk.atmcrds)
						if (rp_c.distance(lk_c) < 1.5)
						{
							clash = true;
							break;
						}
					if (clash) break;
				}
				if (clash) break;
			}
		}
		if (!clash)
		{
			std::string new_name = "new_Pocket_" + to_string(i) + ".pdb";
			pocket.writepocket(new_name);
			mi_new.insert(std::make_pair(itm->first, mi)); // <score, mi>
			i++;
		}
	}
	return mi_new;
}

std::map<double, MatchInfo> Theozyme::double_check(int outputmatchnum, std::map<double, MatchInfo> mis,
		std::map<string, vector<vector<XYZ>>> sca_crds_all,
		std::map<string, vector<KeyAtom>> sca_kas,
		double plmc_clash, double prpr_clash, double ls_clash, vector<int> rotlibc)
{
	ofstream ofs("clash.txt");
	std::map<double, MatchInfo> mis_new;
	int i = 0;
	vector<vector<string>> nr_rcs; // judge if mi is non-redundant. <<cid_rid>>
	int count = 1;
	for (auto itm = mis.begin(); itm != mis.end(); itm++)
	{
		if (i >= outputmatchnum) break;
		if (count % 1000 == 0) cout << "Finish checking " << count << " matches." << endl;
		count++;
		MatchInfo mi = itm->second;
		// sca_mc;
		string sf = mi.scaffold_filename;
		vector<vector<XYZ>> sca_mc = sca_crds_all[sf];
		auto sk = sca_kas[mi.scaffold_filename];

		bool isnew = true;
		vector<string> rcs;
		for (auto si : mi.matchposes)
			if (si != -1)
				rcs.push_back(sk[si].cid + "_" + to_string(sk[si].rid));
			else
				rcs.push_back("Nomatch");
		if (nr_rcs.size() == 0)
			nr_rcs.push_back(rcs);
		else
		{
			for (auto nn : nr_rcs)
			{
				int same = 0;
				for (int i = 0; i < mi.matchposes.size(); i++)
					if (nn[i] == rcs[i]) same++;
				if (same == mi.matchposes.size()) isnew = false;
				if (!isnew) break;
			}
		}
		if (!isnew) continue;

// checking clash
		// poc_lig & mc;
		string pf = mi.pocket_filename;
		Pocket pocket;
		pocket.readrefpdb(pf);
		auto lig_test = pocket.ligands();
		if (!itm->second.dm)
		{
			auto qf = itm->second.qf;
			pocket.transform(qf);
		}
		auto poc_ligs = pocket.ligands();
		bool lig_clash = false;
		for (auto lig : poc_ligs)
		{
			for (auto la : lig.atmcrds)
			{
				for (auto sm : sca_mc)
				{
					for (auto s : sm)
					{
						if (la.distance(s) < plmc_clash)
						{
							ofs << "Clash: Ligand " << la.x_ << " " << la.y_ << " " << la.z_
									<< " with Scaffold_MC " << s.x_ << " " << s.y_ << " " << s.z_
									<< " with RMSD: " << la.distance(s) << endl;
 							lig_clash = true;
						}
						if (lig_clash) break;
					}
					if (lig_clash) break;
				}
				if (lig_clash) break;
			}
			if (lig_clash) break;
		} // lig-ScaMC
		if (!lig_clash)	// lig_mc no clash, then update those rotamers
		{
			// change_sc;
			renew_residues(pocket, sca_kas[mi.scaffold_filename], mi, rotlibc);
			// judge clash between res_sc-res_sc (2.0 A) and lig-res_sc (1.3 A in case that forming covalent bond).
			bool clash = false;
			auto residues = pocket.residues();
			auto ligands = pocket.ligands();
			for (int p = 0; p < residues.size()-1; p++)
			{
				auto rp = residues[p];
				if (mi.matchposes[p] < 0) continue;
				for (int q = p+1; q < residues.size(); q++)
				{
					auto rq = residues[q];
					if (mi.matchposes[q] < 0) continue;
					for (int i = 0; i < rp.atmcrds.size(); i++)
					{
						if (rp.atmnames[i] == "N" || rp.atmnames[i] == "CA" ||
								rp.atmnames[i] == "O" || rp.atmnames[i] == "C")
							continue;
						auto rp_c = rp.atmcrds[i];
						for (int j = 0; j < rq.atmcrds.size(); j++)
						{
							if (rq.atmnames[j] == "N" || rq.atmnames[j] == "CA" ||
									rq.atmnames[j] == "O" || rq.atmnames[j] == "C")
								continue;
							auto rq_c = rq.atmcrds[j];
							if (rp_c.distance(rq_c) < prpr_clash)
							{
								ofs << "Clash: Pocket_SC " << rp_c.x_ << " " << rp_c.y_ << " " << rp_c.z_
										<< " with Pocket_SC " << rq_c.x_ << " " << rq_c.y_ << " " << rq_c.z_
										<< " with RMSD: " << rp_c.distance(rq_c) << endl;
								clash = true;
								break;
							}
						}
						if (clash) break;
					}
					if (clash) break;
				}
				if (clash) break;
			} // PocSC-PocSC
			for (int p = 0; p < residues.size(); p++)
			{
				auto rp = residues[p];
				if (mi.matchposes[p] < 0) continue;
				for (int k = 0; k < ligands.size(); k++)
				{
					auto lk = ligands[k];
					for (int i = 0; i < rp.atmcrds.size(); i++)
					{
						if (rp.atmnames[i] == "N" || rp.atmnames[i] == "CA" ||
								rp.atmnames[i] == "O" || rp.atmnames[i] == "C")
							continue;
						auto rp_c = rp.atmcrds[i];
						for (auto lk_c : lk.atmcrds)
							if (rp_c.distance(lk_c) < ls_clash)
							{
								ofs << "Clash: Pocket_SC " << rp_c.x_ << " " << rp_c.y_ << " " << rp_c.z_
										<< " with Ligand " << lk_c.x_ << " " << lk_c.y_ << " " << lk_c.z_
										<< " with RMSD: " << rp_c.distance(lk_c) << endl;
								clash = true;
								break;
							}
						if (clash) break;
					}
					if (clash) break;
				}
				if (clash) break;
			} // Lig-PocSC
			if (!clash)
			{
				std::string new_name = "new_Pocket_" + to_string(i) + ".pdb";
				pocket.writepocket(new_name);
				mis_new.insert(std::make_pair(itm->first, mi)); // <score, mi>
				i++;
				nr_rcs.push_back(rcs); // new rcs & no clash, update nr_rcs.
			}
		}
	}
	return mis_new;
}

std::map<double, MatchInfo> Theozyme::checkclash_mis(std::map<double, MatchInfo> mis,
		std::map<std::string, std::vector<std::vector<XYZ>>> sca_crds_all, vector<double> RMSD_ths)
{
	std::map<double, MatchInfo> mis_new;
	double r_max = 0.0;
	for (auto r : RMSD_ths)
		if (r > r_max)
			r_max = r;
	for (auto itm = mis.begin(); itm != mis.end(); itm++)
	{
		MatchInfo mi = itm->second;
		// copy mis_member
		auto mp = mi.matchposes;
		auto rmsds = mi.rmsds;

		// sca_mc;
		string sf = mi.scaffold_filename;
		vector<vector<XYZ>> sca_mc = sca_crds_all[sf];

		// poc_sc;
		string pf = mi.pocket_filename;
		Pocket pocket;
		pocket.readrefpdb(pf);
		auto qf = itm->second.qf;
		pocket.transform(qf);
		auto poc_res = pocket.residues();
		for (int i = 0; i < poc_res.size(); i++)
		{
			int clash_num = 0;
			vector<XYZ> poc_sc;
			for (int j = 0; j < poc_res[i].atmnames.size(); j++)
			{
				string an = poc_res[i].atmnames[j];
				if (an != "C" && an != "O" && an != "N") // CA should also be considered for its pos.
					poc_sc.push_back(poc_res[i].atmcrds[j]); // sc_crds;
			}
			for (auto sm : sca_mc) // each sca_res
			{
				bool clash = false;
				for (auto s : sm) // each s_mc_atm
				{
					if (clash) break;
					for (auto p : poc_sc) // each p_sc_atm
						if (s.distance(p) < 1.6) // 2.0 may be to strict
						{
							clash = true;
							break;
						}
				}
				if (clash) clash_num++;
			}
			if (clash_num > 1) // 1 is for its pos.
				rmsds[i] = r_max;
		} // each poc_res

		mi.rmsds = rmsds;
		// mis_new:insert
		double ave_r = 0.0;
		for (int r = 0; r < rmsds.size(); r++)
			ave_r += rmsds[r]/(RMSD_ths[r]*RMSD_ths[r]);
		mis_new.insert(std::make_pair(100*ave_r/rmsds.size(), mi)); // <score, mi>
	} // each match
	return mis_new;
}

AtomList Theozyme::readal(string filename)
{
	std::map<string, int> record_residue;
    std::ifstream ifs(filename);
    if (!ifs.good())
    {
    	cout << "AtomList reading failure." << endl;
    	exit(1);
    }
    AtomList al;
    while (true)
    {
		std::string line;
		line.clear();
		getline(ifs, line);
		if(!ifs.good()) break;
		if (line.size()==0) continue;
		std::vector<std::string> words;
		std::stringstream input(line);
		std::string word;
		while(input>>word) words.push_back(word);
		assert(words.size() > 3);
    	al.cids.push_back(words[0]);
    	al.rids.push_back(words[1]);
    	al.rname.push_back(words[2]);
    	if (!isresidue(words[2]))
    	{
    		cout << "Error in AtomList: residue " << words[2] << " is not a common residue." << endl;
    		exit(1);
    	}
    	vector<string> new_aname;
    	for (int i = 3; i < words.size(); i++)
    		new_aname.push_back(words[i]);
    	al.anames.push_back(new_aname);
    	if (record_residue.count(words[2]) > 0)
    	{
    		record_residue[words[2]] = record_residue[words[2]] + 1;
    		al.ids.push_back(record_residue[words[2]]);
    	}
    	else
    	{
    		record_residue.insert(make_pair(words[2], 1));
    		al.ids.push_back(1);
    	}
    }
    return al;
}

void Theozyme::readrotlib(vector<int> &rotlibc, string filename)
{
    std::ifstream ifs(filename);
    if (!ifs.good())
    {
    	cout << "RotLibChoice reading failure." << endl;
    	exit(1);
    }
    rotlibc.assign(100, 1);
    while (true)
    {
		std::string line;
		line.clear();
		getline(ifs, line);
		if(!ifs.good()) break;
		if (line.size()==0) continue;
		std::vector<std::string> words;
		std::stringstream input(line);
		std::string word;
		while(input>>word) words.push_back(word);
		if (words[0] == "default" && words[2] == "for")
			rotlibc.assign(stoi(words[3]), stoi(words[1]));
		if (words.size() == 2)
		{
			if (rotlibc.size() > stoi(words[0]))
				rotlibc[stoi(words[0])] = stoi(words[1]);
			else
			{
				cout << "Error: " << words[0] << " is larger than the size of rotlibchoice, which is "
						<< rotlibc.size() << endl;
			}
		}
    }
    cout << "UserDefined rotLib:" << endl;
    for (int i = 0; i < rotlibc.size(); i++)
    	cout << " Pocket Residue " << i << " using rotLib mode " << rotlibc[i] << endl;
}

vector<KeyAtom> Theozyme::extractka(AtomList al, Pocket poc)
{
	vector<KeyAtom> kas;
	auto residues = poc.residues();
	if (residues.size() != al.rname.size())
	{
		cout << "PocketResidues (" << residues.size() << ") are not correctly described in al.txt (" << al.rname.size() << "). (in residue number)" << endl;
		exit(1);
	}
	int seq_old = -1; // used to check if al.txt is written in residue sequence.
	for (int i = 0; i < al.cids.size(); i++)
	{
		int seq_new = 0;
		bool find_res = false;
		for (auto res : residues)
		{
			seq_new++;
			if (res.chainid == al.cids[i] && to_string(res.resid) == al.rids[i])
			{
				if (seq_new < seq_old)
				{
					cout << "AtomList file should be written in the same sequence as pocket residues." << endl;
					exit(1);
				}
				else
					seq_old = seq_new;
				KeyAtom ka;
				ka.rname = res.resname;
				ka.id = al.ids[i];
				ka.rid = res.resid;
				ka.cid = res.chainid;
				for (auto a : al.anames[i])
				{
					bool find_atm = false;
					for (int j = 0; j < res.atmnames.size(); j++)
						if (a == res.atmnames[j])
						{
							ka.crds.push_back(res.atmcrds[j]);
							ka.anames.push_back(a);
							find_atm = true;
							break;
						}
					if (!find_atm)
					{
						cout << "AtomName " << a << " is not in the residue with chainID " << al.cids[i]
					        << " resID " << al.rids[i] << endl;
						exit(1);
					}
				} // pocket_atm
				find_res = true;
				kas.push_back(ka);
			}
			if (find_res) break;
		} // pocket_res.
		if (!find_res)
		{
			cout << "No such residue with chainID " << al.cids[i]
		        << " resID " << al.rids[i] << " in the Pocket." << endl;
			exit(1);
		}
	}
	return kas;
}

void Theozyme::readkr_rmsdths(vector<double> &RMSD_ths, string filename)
{
	ifstream ifk;
	ifk.open(filename.c_str());
	if(!ifk.good()) {
		std::cout << "kr.txt file failure" << std::endl;
		exit(1);
	}
	double kr_th = 2.0; // default RMSD_th for key residues.
	while(true) {
		std::string line;
		line.clear();
		getline(ifk, line);
		if(!ifk.good()) break;
		if (line.size()==0) continue;
		std::vector<std::string> words;
		std::stringstream input(line);
		std::string word;
		while(input>>word) words.push_back(word);
		if (words.size() == 2 && words[0] == "RMSD")
			kr_th = std::stod(words[1]);
		else if (words.size() == 1)
		{
			assert(RMSD_ths.size() > std::stoi(words[0])-1);
			RMSD_ths[std::stoi(words[0])-1] = kr_th;
		}
		else
		{
			std::cout << "kr.txt should only contain RMSD value & residue_id";
			exit(1);
		}
	}
	ifk.close();
}

void Theozyme::readkr_rmsdths(vector<double> &RMSD_ths, string filename, set<int> &kr_ids)
{
	ifstream ifk;
	ifk.open(filename.c_str());
	if(!ifk.good()) {
		std::cout << "kr.txt file failure" << std::endl;
		exit(1);
	}
	double kr_th = 2.0; // default RMSD_th for key residues.
	while(true) {
		std::string line;
		line.clear();
		getline(ifk, line);
		if(!ifk.good()) break;
		if (line.size()==0) continue;
		std::vector<std::string> words;
		std::stringstream input(line);
		std::string word;
		while(input>>word) words.push_back(word);
		if (words.size() == 2 && words[0] == "RMSD")
			kr_th = std::stod(words[1]);
		else if (words.size() == 1)
		{
			if(RMSD_ths.size() <= std::stoi(words[0]))
			{
				cout << "Key Residue file contains Residue_id exceeding the residue number in the pocket." << endl;
				exit(1);
			}
			RMSD_ths[std::stoi(words[0])] = kr_th;
			kr_ids.insert(std::stoi(words[0]));
		}
		else
		{
			std::cout << "kr.txt should only contain RMSD value & residue_id";
			exit(1);
		}
	}
	ifk.close();
}

pair<float, float> Theozyme::getppfromgrs(GeneralRes pgr, GeneralRes gr, GeneralRes ngr)
{
	pair<float, float> pp;

	// calculate phi & psi
	XYZ ncrd, ccrd, cacrd, pccrd, nncrd;
	for (int a = 0; a < gr.atmcrds.size(); a++)
		if (gr.atmnames[a] == "N") ncrd = gr.atmcrds[a];
		else if (gr.atmnames[a] == "CA") cacrd = gr.atmcrds[a];
		else if (gr.atmnames[a] == "C") ccrd = gr.atmcrds[a];
	for (int a = 0; a < pgr.atmcrds.size(); a++)
		if (pgr.atmnames[a] == "C") pccrd = pgr.atmcrds[a];
	for (int a = 0; a < ngr.atmcrds.size(); a++)
		if (ngr.atmnames[a] == "N") ncrd = ngr.atmcrds[a];
	float rad = 180.0 / 3.14159265358979323846;
	float phi = NSPgeometry::torsion(pccrd, ncrd, cacrd, ccrd);
    phi *= rad;
    while (phi > 180.0)
        phi -= 360.0;
    while (phi < -180.0)
        phi += 360.0;
    float psi = NSPgeometry::torsion(ncrd, cacrd, ccrd, nncrd);
    psi *= rad;
    while (psi > 180.0)
        psi -= 360.0;
    while (psi < -180.0)
        psi += 360.0;
    pp.first = phi;
    pp.second = psi;
	return pp;
}

vector<vector<XYZ>> Theozyme::get_mc_crds(vector<GeneralRes> residues)
{
	vector<vector<XYZ>> mc_crds;
	for (auto r : residues)
	{
		vector<XYZ> crd;
		for (int i = 0; i < r.atmnames.size(); i++)
			if (r.atmnames[i] == "N" || r.atmnames[i] == "CA" ||
					r.atmnames[i] == "C" || r.atmnames[i] == "O")
				crd.push_back(r.atmcrds[i]);
		mc_crds.push_back(crd);
	}
	return mc_crds;
}

bool Theozyme::mc_rg_hasclash(vector<XYZ> sc_crds, vector<vector<XYZ>> sca_mc_crds, int thissite)
{
	bool hasclash = false;
	for (int i = 0; i < sca_mc_crds.size(); i++)
	{
		if (i == thissite)
			continue;
		for (auto sc : sc_crds)
		{
			for (auto mc : sca_mc_crds[i])
			{
				double d = sc.distance(mc);
				if (d > 15)
					continue;
				else if (d < 2.5)
					hasclash = true;
				if (hasclash) break;
			}
			if (hasclash) break;
		}
		if (hasclash) break;
	}
	return hasclash;
}

void Theozyme::renew_residues(Pocket &poc, vector<KeyAtom> sca_ka, MatchInfo mi, vector<int> rotlibc)
{
	auto residues = poc.residues();
	for (int gid = 0; gid < residues.size(); gid++)
	{
		if (mi.matchposes[gid] < 0 || mi.matchposes[gid] >= sca_ka.size()) // not match or wrong
			continue;
		KeyAtom sk = sca_ka[mi.matchposes[gid]];
		std::map<string, XYZ> aname_crd;
		aname_crd.insert(make_pair("N", sk.ncaco[0]));
		aname_crd.insert(make_pair("CA", sk.ncaco[1]));
		aname_crd.insert(make_pair("C", sk.ncaco[2]));
		aname_crd.insert(make_pair("O", sk.ncaco[3]));

		PhipsiLib ppLib;
		XYZ thisncrd, thisccrd, thiscacrd, can, cac, x, y, z;
		thisncrd = aname_crd["N"];
		thisccrd = aname_crd["C"];
		thiscacrd = aname_crd["CA"];
		can = ~(thisncrd-thiscacrd);
		cac = ~(thisccrd-thiscacrd);
		x = ~(can+cac);
		z = ~(can^cac);
		y = ~(z^x);
		LocalFrame cs;
		cs.origin_ = thiscacrd;
		cs.axis_.push_back(x);
		cs.axis_.push_back(y);
		cs.axis_.push_back(z);
		Phipsi pp(sk.phi, sk.psi);
		int ppType = ppLib.phipsiToIndex(&pp);
		RotamerLib* rotLib;
		if (rotlibc.size() > 0)
		{
			switch(rotlibc[gid])
			{
			case 0:
			{
				if(pp.regionAB() == 'A') rotLib = new RotamerLib("A0");
				else rotLib = new RotamerLib("B0");
				break;
			}
			case 1:
			{
				if(pp.regionAB() == 'A') rotLib = new RotamerLib("A1");
				else rotLib = new RotamerLib("B1");
				break;
			}
			case 2:
			{
				if(pp.regionAB() == 'A') rotLib = new RotamerLib("A2");
				else rotLib = new RotamerLib("B2");
				break;
			}
			case 3:
			{
				if(pp.regionAB() == 'A') rotLib = new RotamerLib("A3");
				else rotLib = new RotamerLib("B3");
				break;
			}
			}
		}
		else
		{
			if(pp.regionAB() == 'A') rotLib = new RotamerLib("A1");
			else rotLib = new RotamerLib("B1");
		}
		RotamerGroup* gp = rotLib->getAAGroup(sk.rname); // rotlists for certain res.
		for (auto rl : gp->rotList)
		{
			// judge if it is the wanted rotamer.
			vector<XYZ> scTerms;
			rl->buildSidechain(cs, scTerms);
			bool find = true;
			int f = 0;
			for (int s = 0; s < sk.crds.size(); s++)
			{
				if (sk.anames[s] != "N" && sk.anames[s] != "CA" && sk.anames[s] != "C" && sk.anames[s] != "O")
				{
					f++;
					auto kacrd = sk.crds[s];
					for (auto sT : scTerms)
						if (sT.distance(kacrd) < 0.01)
						{
							f--;
							break;
						}
				} // only check SC_atoms
			}
			if (f != 0) find = false;
			if (find)
			{
				Residue* res = new Residue();
				res->addAtom(new Atom("N", thisncrd));
				res->addAtom(new Atom("CA", thiscacrd));
				res->addAtom(new Atom("C", thisccrd));
				res->addAtom(new Atom("O", sk.ncaco[3]));
				res->buildRotamer(rl);
	    		for (int i = 0; i < res->getAtomList()->size(); i++)
	    		{
	    			 auto atm = res->getAtomList()->at(i);
	    			 aname_crd.insert(std::make_pair(atm->name, atm->getCoord()));
	    		}
	    		// update Pocket
	    		auto residue = residues[gid];
	    		vector<XYZ> newcrds;
	    		for (auto ra : residue.atmnames)
	    			newcrds.push_back(aname_crd[ra]);
	    		poc.crdchange(false, gid, newcrds);
				break;
			}
		} // search for that rotamer
	}
}

vector<KeyAtom> Theozyme::designka(double erotcutoff, Pocket sca, AtomList al, Region reg, vector<int> rotlibc)
{
	vector<KeyAtom> ka;
	PhipsiLib ppLib;
	auto residues = sca.residues();
	vector<vector<XYZ>> mc_crds = get_mc_crds(residues);
	for (int r = 1; r < residues.size()-1; r++) // first and last residue will be skipped
	{
		auto res = residues[r];
		for (int i = 0; i < reg.cids.size(); i++)
		{
			if (res.chainid == reg.cids[i] &&
					res.resid >= reg.i_begins[i] &&
					res.resid <= reg.i_ends[i])
			{
				XYZ thisncrd, thisccrd, thiscacrd, can, cac, x, y, z, thisocrd;
				for (int a = 0; a < res.atmnames.size(); a++)
					if (res.atmnames[a] == "N") thisncrd = res.atmcrds[a];
					else if (res.atmnames[a] == "CA") thiscacrd = res.atmcrds[a];
					else if (res.atmnames[a] == "C") thisccrd = res.atmcrds[a];
					else if (res.atmnames[a] == "O") thisocrd = res.atmcrds[a];
				can = ~(thisncrd-thiscacrd);
				cac = ~(thisccrd-thiscacrd);
				x = ~(can+cac);
				z = ~(can^cac);
				y = ~(z^x);
				LocalFrame cs;
				cs.origin_ = thiscacrd;
				cs.axis_.push_back(x);
				cs.axis_.push_back(y);
				cs.axis_.push_back(z);
				pair<float, float> pnp = getppfromgrs(residues[r-1], res, residues[r+1]);
				Phipsi pp(pnp.first, pnp.second);
				int ppType = ppLib.phipsiToIndex(&pp);
				RotamerLib* rotLib;
//				if(pp.regionAB() == 'A') rotLib = new RotamerLib("A1");
//				else rotLib = new RotamerLib("B1");
				// collect gp
				for (int a = 0; a < al.rname.size(); a++)
				{
					if (rotlibc.size() > 0)
					{
						switch(rotlibc[a])
						{
						case 0:
						{
							if(pp.regionAB() == 'A') rotLib = new RotamerLib("A0");
							else rotLib = new RotamerLib("B0");
							break;
						}
						case 1:
						{
							if(pp.regionAB() == 'A') rotLib = new RotamerLib("A1");
							else rotLib = new RotamerLib("B1");
							break;
						}
						case 2:
						{
							if(pp.regionAB() == 'A') rotLib = new RotamerLib("A2");
							else rotLib = new RotamerLib("B2");
							break;
						}
						case 3:
						{
							if(pp.regionAB() == 'A') rotLib = new RotamerLib("A3");
							else rotLib = new RotamerLib("B3");
							break;
						}
						}
					}
					else // default: A1 or B1
					{
						if(pp.regionAB() == 'A') rotLib = new RotamerLib("A1");
						else rotLib = new RotamerLib("B1");
					}
					RotamerGroup* gp = rotLib->getAAGroup(al.rname[a]);
					for (int i = 0; i < gp->rotNum; i++)
					{
						Rotamer* rot = gp->rotList[i];
						float eRot = rotLib->getRotamerEnergy(rot->rotName, ppType);
						if(eRot >= erotcutoff) // delete high_Erot: 4 (ABACUS). here we use 5
						{
							gp->deleteRotamer(i);
							i--;
						}
						else // delete sc_mc clash
						{
							vector<XYZ> scTerms;
							rot->buildSidechain(cs, scTerms);
							bool hasclash = mc_rg_hasclash(scTerms, mc_crds, r);
							if (hasclash)
							{
								gp->deleteRotamer(i);
								i--;
							}
						}
					}
					// gp->ka
					for (auto aarot : gp->aaRots)
					{
						for (auto aar : *aarot)
						{
							KeyAtom k;
							k.rname = al.rname[a];
							k.id = al.ids[a];
							k.cid = res.chainid;
							k.rid = res.resid;
							for (auto an : al.anames[a])
							{
								if (an == "N")
								{
									k.anames.push_back(an);
									k.crds.push_back(thisncrd);
								}
								else if (an == "CA")
								{
									k.anames.push_back(an);
									k.crds.push_back(thiscacrd);
								}
								else if (an == "C")
								{
									k.anames.push_back(an);
									k.crds.push_back(thisccrd);
								}
								else if (an == "O")
								{
									k.anames.push_back(an);
									k.crds.push_back(thisocrd);
								}
								else
								{
									bool find_atm = false;
									for (int c = 0; c < aar->coordList.size(); c++)
										if (an == aar->atomNameList[c])
										{
											k.anames.push_back(an);
											k.crds.push_back(cs.local2globalcrd(aar->coordList[c]));
											find_atm = true;
											break;
										} // sc_atom
									if (!find_atm)
									{
										cout << "Couldn't find " << an << " in the created " << al.rname[a] <<
												" in the scaffold." << endl;
										exit(1);
									}
								}
							}
							k.phi = pp.phi;
							k.psi = pp.psi;
							k.ncaco.push_back(thisncrd);
							k.ncaco.push_back(thiscacrd);
							k.ncaco.push_back(thisccrd);
							k.ncaco.push_back(thisocrd);
							ka.push_back(k);
						} // each possible rot
					}
				} // each kr_type
			} // each region site.
		} // cid of region.
	}

	return ka;
}

KeyAtom Theozyme::extractka_fromrot(AtomList al, Rotamer* rot)
{
	KeyAtom k;
	for (int i = 0; i < al.rname.size(); i++)
	{
		if (al.rname[i] == rot->triName)
		{
			k.rname = rot->triName;
			for (int a = 0; a < rot->atomNameList.size(); a++)
			{
				k.anames.push_back(rot->atomNameList[a]);
				k.crds.push_back(rot->coordList[a]);
				k.id = al.ids[i];
			}
		}
	}

	return k;
}

char Theozyme::tri2one(string tri)
{
	std::map<string, char> t2o;
	t2o.insert(make_pair("GLY", 'G'));
	t2o.insert(make_pair("ALA", 'A'));
	t2o.insert(make_pair("VAL", 'V'));
	t2o.insert(make_pair("LEU", 'L'));
	t2o.insert(make_pair("ILE", 'I'));
	t2o.insert(make_pair("MET", 'M'));
	t2o.insert(make_pair("TRP", 'W'));
	t2o.insert(make_pair("PHE", 'F'));
	t2o.insert(make_pair("PRO", 'P'));
	t2o.insert(make_pair("SER", 'S'));
	t2o.insert(make_pair("THR", 'T'));
	t2o.insert(make_pair("CYS", 'C'));
	t2o.insert(make_pair("TYR", 'Y'));
	t2o.insert(make_pair("ASN", 'N'));
	t2o.insert(make_pair("GLN", 'Q'));
	t2o.insert(make_pair("ASP", 'D'));
	t2o.insert(make_pair("GLU", 'E'));
	t2o.insert(make_pair("LYS", 'K'));
	t2o.insert(make_pair("ARG", 'R'));
	t2o.insert(make_pair("HIS", 'H'));
	return t2o[tri];
}

Region Theozyme::readregion(string region)
{
	Region reg;
	ifstream ifs(region);
	if(!ifs.good()) {
		std::cout << "region file failure" << std::endl;
		exit(1);
	}
	while(true) {
		std::string line;
		line.clear();
		getline(ifs, line);
		if(!ifs.good()) break;
		if (line.size()==0) continue;
		std::vector<std::string> words;
		std::stringstream input(line);
		std::string word;
		while(input>>word) words.push_back(word);
		reg.cids.push_back(words[0]);
		reg.i_begins.push_back(stoi(words[1]));
		reg.i_ends.push_back(stoi(words[2]));
	}
	ifs.close();
	return reg;
}

bool Theozyme::passscreen(Pocket poc, int cn, int maxn, int minn)
{
	int chainnum = 1;
	int aanum = 0;
	bool passmin = false; // at least one chain should have >=minn AANum.
	auto res = poc.residues();
	if (res.size() == 0) return false; // should have at least one protein chain.
	string oldcid = res[0].chainid;
	for (auto r : res)
		if (r.chainid != oldcid)
		{
			chainnum++;
			oldcid = r.chainid;
			aanum = 1;
		}
		else
		{
			aanum++;
			if (aanum > maxn) return false;
			if (!passmin && aanum >= minn) passmin = true;
		}
	if (!passmin) return false;
	if (chainnum > cn) return false;
	else return true;
}

bool Theozyme::passscreen(Pocket poc, int cn, int maxn, int minn, bool NoMetal, bool NoApo, bool NoBig)
{
	int chainnum = 1;
	int aanum = 0;
	bool passmin = false; // at least one chain should have >=minn AANum.
	auto res = poc.residues();
	if (res.size() == 0) return false; // should have at least one protein chain.
	string oldcid = res[0].chainid;
	for (auto r : res)
		if (r.chainid != oldcid)
		{
			chainnum++;
			oldcid = r.chainid;
			aanum = 1;
		}
		else
		{
			aanum++;
			if (aanum > maxn) return false;
			if (!passmin && aanum >= minn) passmin = true;
		}
	if (!passmin) return false;
	if (chainnum > cn) return false;

	// NoMetal NoApo NoBig
	auto lig = poc.ligands();
	if (NoMetal || NoApo || NoBig)
	{
		bool havinglig = false;
		for (auto l : lig)
		{
			if (l.resname.size() == 2) // could be metal or DNA
			{
				if (NoMetal && l.atmnames.size() == 1)
					return false;
				if (NoBig)
					if (l.resname == "DA" || l.resname == "DT" ||
							l.resname == "DC" || l.resname == "DG")
						return false;
			}
			else if (l.resname.size() == 1 && NoBig) // RNA
				return false;
			else if (l.resname.size() == 3) // HOH, PX4, ligand.
			{
				if (l.resname != "HOH" && l.resname != "PX4" &&
						l.resname != "MSE" && l.resname != "SEC")
					havinglig = true;
				if (l.resname == "PX4" && NoBig)
					return false;
			}
		}
		if (NoApo && !havinglig)
			return false;
	}

	return true;
}

map<int, vector<int>> Theozyme::rankligmw(vector<Pocket> pocs, vector<string> pfns, int lmw)
{
	map<int, vector<int>> rankmaps;
	for (int i = 0; i < pfns.size(); i++)
	{
		int diff = 9999;
		auto ligs = pocs[i].ligands();
		for (auto l : ligs)
		{
			int lmw_new = 0;
			for (auto an : l.atmnames)
				if (an[0] != 'H') lmw_new++;
			diff = abs(lmw_new - lmw) < diff ? abs(lmw_new - lmw) : diff;
		}
		if (rankmaps.count(diff) > 0)
			rankmaps[diff].push_back(i);
		else
		{
			vector<int> grp{i};
			rankmaps.insert(make_pair(diff, grp));
		}
	}
	return rankmaps;
}

map<int, vector<int>> Theozyme::rankligmw(vector<Pocket> pocs, vector<string> pfns, int lmw,
		int InnerCA, double Radius)
{
	map<int, vector<int>> rankmaps;
	for (int i = 0; i < pfns.size(); i++)
	{
		int diff = 9999;
		auto ligs = pocs[i].ligands();
		int similig = -1;
		for (int a = 0; a < ligs.size(); a++)
		{
			auto l = ligs[a];
			int lmw_new = 0;
			for (auto an : l.atmnames)
				if (an[0] != 'H') lmw_new++;
			if (abs(lmw_new - lmw) < diff)
			{
				diff = abs(lmw_new - lmw);
				similig = a;
			}
		}
		auto res = pocs[i].residues();
		int canum = 0;
		for (auto r : res)
			for (int ra = 0; ra < r.atmnames.size(); ra++)
			{
				if (r.atmnames[ra] == "CA")
					for (int b = 0; b < ligs[similig].atmcrds.size(); b++)
					{
						auto lcrd = ligs[similig].atmcrds[b];
						if (ligs[similig].atmnames[b][0] == 'H') continue;
						if (lcrd.distance(r.atmcrds[ra]) < Radius)
						{
							canum++;
							break;
						}
					}
			}
		if (canum < InnerCA) diff = diff + 100;

		if (rankmaps.count(diff) > 0)
			rankmaps[diff].push_back(i);
		else
		{
			vector<int> grp{i};
			rankmaps.insert(make_pair(diff, grp));
		}
	}
	return rankmaps;
}

int Theozyme::cubenumber(ScaPocket sp, double cubesize)
{
	std::vector<XYZ> mccrds = sp.pocmcs();
	XYZ center = sp.center();
	double radi = sp.radi();
	int step = ceil(radi/cubesize);
	int cubenum = 0;  // cavity should be continuously, so in the following we always count from center.
	for (int x = 0; x <= step; x++)
	{
		for (int y = 0; y <= step; y++)
		{
			for (int z = 0; z <= step; z++)
			{
				XYZ cubecenter {center.x_+x*cubesize, center.y_+y*cubesize, center.z_+z*cubesize};

			}
		}
	}
}

//void ScaPocket::poc2scapocket(string filename, Pocket poc, string ligname, double radi)
//{
//
//}

void ScaPocket::poc2scapocket(string filename, Pocket poc, XYZ center, double radi)
{
	name_ = filename;
	center_.x_ = center.x_;
	center_.y_ = center.y_;
	center_.z_ = center.z_;
	radi_ = radi;
	auto res = poc.residues();
	for (auto r : res)
		for (int i = 0; i < r.atmnames.size(); i++)
		{
			if (r.atmnames[i] == "N" || r.atmnames[i] == "CA" ||
					r.atmnames[i] == "C" || r.atmnames[i] == "O")
				if (r.atmcrds[i].distance(center) < radi)
					pocmcs_.push_back(r.atmcrds[i]);
		}
}

bool Theozyme::withindistance(GeneralRes a, GeneralRes b, double d)
{
	for (auto acrd : a.atmcrds)
		for (auto bcrd : b.atmcrds)
			if (acrd.distance(bcrd) < d)
				return true;
	return false;
}

map<string, vector<int>> Theozyme::findclashRotLib(double erotcutoff, vector<GeneralRes> ms_rs, vector<GeneralRes> mp_rs,
		vector<GeneralRes> mp_li, map<int, int> mpsites, double pc_th, double npc_th)
{
	map<string, vector<int>> clashsites_ms;
	PhipsiLib ppLib;
	vector<vector<XYZ>> mc_crds = get_mc_crds(ms_rs);
	set<int> sites;
	for (auto iter = mpsites.begin(); iter != mpsites.end(); iter++)
		sites.insert(iter->second);
	for (int r = 1; r < ms_rs.size()-1; r++) // first and last residue will be skipped
	{
		if (sites.count(r) > 0) continue; // matched sites don't need checkclash.
		auto res = ms_rs[r];
		if (res.resname == "GLY") continue; // GLY shouldn't have clash.
		// if Ca_ms_rs - pocket < 10.0, then do the following.
		bool within = false;
		for (int p = 0; p < mp_rs.size(); p++)
		{
			if (mpsites.count(p) > 0)
				if (withindistance(res, mp_rs[p], 4.5)) within = true;
		} // no clash wiht matched PocRes.
		for (auto mpli : mp_li)
			if (withindistance(res, mpli, 4.5)) within = true;
		if (within)
		{
			XYZ thisncrd, thisccrd, thiscacrd, can, cac, x, y, z, thisocrd;
			for (int a = 0; a < res.atmnames.size(); a++)
				if (res.atmnames[a] == "N") thisncrd = res.atmcrds[a];
				else if (res.atmnames[a] == "CA") thiscacrd = res.atmcrds[a];
				else if (res.atmnames[a] == "C") thisccrd = res.atmcrds[a];
				else if (res.atmnames[a] == "O") thisocrd = res.atmcrds[a];
			can = ~(thisncrd-thiscacrd);
			cac = ~(thisccrd-thiscacrd);
			x = ~(can+cac);
			z = ~(can^cac);
			y = ~(z^x);
			LocalFrame cs;
			cs.origin_ = thiscacrd;
			cs.axis_.push_back(x);
			cs.axis_.push_back(y);
			cs.axis_.push_back(z);
			pair<float, float> pnp = getppfromgrs(ms_rs[r-1], res, ms_rs[r+1]);
			Phipsi pp(pnp.first, pnp.second);
			int ppType = ppLib.phipsiToIndex(&pp);
			RotamerLib* rotLib;
			if(pp.regionAB() == 'A') rotLib = new RotamerLib("A1");
			else rotLib = new RotamerLib("B1");
			RotamerGroup* gp = rotLib->getAAGroup(res.resname);
			bool clash = true;
			for (int i = 0; i < gp->rotNum; i++)
			{
				if (!clash) break;
				Rotamer* rot = gp->rotList[i];
				float eRot = rotLib->getRotamerEnergy(rot->rotName, ppType);
				if(eRot < erotcutoff)
				{
					vector<XYZ> scTerms;
					rot->buildSidechain(cs, scTerms);
//					bool hasclash = mc_rg_hasclash(scTerms, mc_crds, r);
					// check if clash with pocket.
					bool findclash = false;
					for (int s = 0; s < scTerms.size(); s++)
					{
						if (findclash) break;
						auto sct = scTerms[s];
						for (int p = 0; p < mp_rs.size(); p++)
						{
							if (mpsites.count(p) > 0)
							{
								for (int a = 0; a < mp_rs[p].atmcrds.size(); a++)
								{
									if (mp_rs[p].atmnames[a] == "N" || mp_rs[p].atmnames[a] == "CA" ||
											mp_rs[p].atmnames[a] == "C" || mp_rs[p].atmnames[a] == "O")
										continue;
									auto rscrd = mp_rs[p].atmcrds[a];
									if (sct.distance(rscrd) < pc_th)
										findclash = true;
									else if (sct.distance(rscrd) < npc_th)
									{
										if (mp_rs[p].elementnames[a] == "C" && rot->atomNameList[s][0] == 'C')
											findclash = true;
										if (mp_rs[p].elementnames[a] == "C" && rot->atomNameList[s][0] == 'S')
											findclash = true;
										if (mp_rs[p].elementnames[a] == "S" && rot->atomNameList[s][0] == 'C')
											findclash = true;
										if (mp_rs[p].elementnames[a] == "S" && rot->atomNameList[s][0] == 'S')
											findclash = true;
									}
								}
							}
						}
						for (auto mpli : mp_li)
							for (int a = 0; a < mpli.atmcrds.size(); a++)
							{
								auto licrd = mpli.atmcrds[a];
								if (sct.distance(licrd) < pc_th)
									findclash = true;
								else if (sct.distance(licrd) < npc_th)
								{
									if (mpli.elementnames[a] == "C" && rot->atomNameList[s][0] == 'C')
										findclash = true;
									if (mpli.elementnames[a] == "C" && rot->atomNameList[s][0] == 'S')
										findclash = true;
									if (mpli.elementnames[a] == "S" && rot->atomNameList[s][0] == 'C')
										findclash = true;
									if (mpli.elementnames[a] == "S" && rot->atomNameList[s][0] == 'S')
										findclash = true;
								}
							}
					} // check clash for each SCAtom
					if (!findclash) clash = false;
				}
			}
			if (clash)
			{
				if (clashsites_ms.count(res.chainid) > 0)
					clashsites_ms[res.chainid].push_back(res.resid);
				else
					clashsites_ms.insert(make_pair(res.chainid, vector<int> {res.resid}));
			}
		}
	}


	return clashsites_ms;
}

bool Theozyme::checkclash_ResSC(GeneralRes a, GeneralRes b, double pc_th, double npc_th)
{
	for (int i = 0; i < a.atmcrds.size(); i++)
		for (int j = 0; j < b.atmcrds.size(); j++)
		{
			if (a.atmnames[i] == "N" || a.atmnames[i] == "CA" || a.atmnames[i] == "C" || a.atmnames[i] == "O")
				continue;
			if (b.atmnames[j] == "N" || b.atmnames[j] == "CA" || b.atmnames[j] == "C" || b.atmnames[j] == "O")
				continue;
			double r = a.atmcrds[i].distance(b.atmcrds[j]);
			if (r < pc_th) return true;
			else if (r < npc_th)
			{
				if (a.elementnames[i] == "C" && b.elementnames[j] == "C") return true;
				if (a.elementnames[i] == "C" && b.elementnames[j] == "S") return true;
				if (a.elementnames[i] == "S" && b.elementnames[j] == "C") return true;
				if (a.elementnames[i] == "S" && b.elementnames[j] == "S") return true;
			}
		}
	return false;
}
