/*
 * extendkeyresidues.cpp
 *
 *  Created on: 2019年8月20日
 *      Author: yxchen
 */
// this fundtion is used to add glys in both ends of key residues
// so that phipsi of key residues will be meaningful.

#include "proteinrep/pdbreader.h"
#include "dstl/randomengine.h"
#include "pdbstatistics/phipsidistr.h"
#include "proteinrep/idealgeometries.h"
#include "geometry/xyz.h"
#include "geometry/calculators.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
using namespace NSPproteinrep;
using namespace NSPgeometry;
double degree = 3.14159265/180.0;

std::vector<XYZ> mccrds(std::vector<PdbRecord> records)
{
	std::vector<XYZ> crds(4);
	for (auto &r : records)
	{
		if (r.atomname == "N")
			crds[0] = NSPgeometry::XYZ{r.x, r.y, r.z};
		if (r.atomname == "CA")
			crds[1] = NSPgeometry::XYZ{r.x, r.y, r.z};
		if (r.atomname == "C")
			crds[2] = NSPgeometry::XYZ{r.x, r.y, r.z};
		if (r.atomname == "O")
			crds[3] = NSPgeometry::XYZ{r.x, r.y, r.z};
	}
	return crds;
}

void genpregly(std::vector<XYZ> &crd_pre, const std::vector<XYZ> &crd_beg, double b_phi, double psi) {
	double omiga = 180.0;
	if (psi > 180.0)
		psi -= 360.0;
	if (psi < -180.0)
		psi += 360.0;
	IdealGeometries & igdat = IdealGeometries::getGlobalInstance();
	crd_pre.resize(4);
	crd_pre[2] = NSPgeometry::InternaltoXYZ(crd_beg[0], crd_beg[1],
			crd_beg[2], igdat.idealLength("C", "N"),
			igdat.idealAngle("C", "N", "CA"), b_phi * degree);
	crd_pre[3] = NSPgeometry::InternaltoXYZ(crd_pre[2], crd_beg[0],
			crd_beg[1], igdat.idealLength("C", "O"),
			igdat.idealAngle("O", "C", "N"), (180.0-omiga) * degree);
	crd_pre[1] = NSPgeometry::InternaltoXYZ(crd_pre[2], crd_beg[0], crd_beg[1],
					igdat.idealLength("CA", "C"),
					igdat.idealAngle("CA", "C", "N"), omiga * degree);
	crd_pre[0] = NSPgeometry::InternaltoXYZ(crd_pre[1], crd_pre[2], crd_beg[0],
				igdat.idealLength("N", "CA"),
				igdat.idealAngle("N", "CA", "C"), psi * degree);
}

void genbackgly(std::vector<XYZ> &crd_bac, const std::vector<XYZ> &crd_end, double e_psi, double phi) {
	double omiga = 180.0;
	if (phi > 180.0)
		phi -= 360.0;
	if (phi < -180.0)
		phi += 360.0;
	IdealGeometries & igdat = IdealGeometries::getGlobalInstance();
	crd_bac.resize(4);
	crd_bac[0] = NSPgeometry::InternaltoXYZ(crd_end[2], crd_end[1],
			crd_end[0], igdat.idealLength("C", "N"),
			igdat.idealAngle("CA", "C", "N"), e_psi * degree);
	crd_bac[1] = NSPgeometry::InternaltoXYZ(crd_bac[0], crd_end[2],
					crd_end[1], igdat.idealLength("N", "CA"),
					igdat.idealAngle("C", "N", "CA"), omiga * degree);
	crd_bac[2] = NSPgeometry::InternaltoXYZ(crd_bac[1], crd_bac[0], crd_end[2],
					igdat.idealLength("CA", "C"),
					igdat.idealAngle("N", "CA", "C"), phi * degree);
	crd_bac[3] = NSPgeometry::InternaltoXYZ(crd_bac[2], crd_bac[1], crd_bac[0],
					igdat.idealLength("C", "O"),
					igdat.idealAngle("CA", "C", "O"), (e_psi + 180.0) * degree);
}

std::vector<PdbRecord> make_pregly(std::vector<XYZ> crds, char chainid)
{
	std::vector<PdbRecord> pregly(4);
	for (int i = 0; i < 4; i++)
	{
		auto &pg = pregly[i];
		if (i == 0) pg.namesymbol = "N";
		if (i == 1 || i == 2) pg.namesymbol = "C";
		if (i == 3) pg.namesymbol = "O";
		if (i == 1) pg.namemodifier = "A";
		pg.label = "ATOM";
		pg.residueid = 0;
		pg.chainid = chainid;
		pg.residuename = "GLY";
		pg.elementname[1] = pg.namesymbol[0];
		pg.x = crds[i].x_;
		pg.y = crds[i].y_;
		pg.z = crds[i].z_;
	}
	return pregly;
}

std::vector<PdbRecord> make_bacgly(std::vector<XYZ> crds, char chainid, int resid)
{
	std::vector<PdbRecord> bacgly(4);
	for (int i = 0; i < 4; i++)
	{
		auto &bg = bacgly[i];
		if (i == 0) bg.namesymbol = "N";
		if (i == 1 || i == 2) bg.namesymbol = "C";
		if (i == 3) bg.namesymbol = "O";
		if (i == 1) bg.namemodifier = "A";
		bg.label = "ATOM";
		bg.residueid = resid;
		bg.chainid = chainid;
		bg.residuename = "GLY";
		bg.elementname[1] = bg.namesymbol[0];
		bg.x = crds[i].x_;
		bg.y = crds[i].y_;
		bg.z = crds[i].z_;
	}
	return bacgly;
}

int main(int argc, char **argv)
{
	// read raw pocket and change to right for
	int n_aid = 1;
	PdbReader pr;
	pr.readpdb(argv[1]);
	std::string chainids = pr.chainids();
	// old_records[chain][res][atm], do not need to record ligand
	std::vector<std::vector<std::vector<PdbRecord>>> o_records(chainids.size()-1);
	// new_records: final form
	std::vector<PdbRecord> n_records;
	std::ofstream ofs("outfile.pdb");
	for (int c = 0; c < chainids.size(); c++)
	{
		char o_cid = chainids[c]; // original chainid
		std::vector<std::string> seq = pr.getaminoacidsequence(o_cid);
		for (int r = 0; r < seq.size(); ++r)
		{
			typename PdbReader::ResKeyType reskey = pr.mappdbkeyint()->pdbResKey(r, c);
			std::vector<PdbRecord> &records = pr.records().at(o_cid).at(reskey);
			if (records[0].label == "HETATM") // ligand
			{
				std::map<std::string, int> atmids; // rename atom names
				for (auto &record : records)
				{
					if (record.label == "HETATM")
					{
						if (record.residuename == "UNL")
						{
							std::string atm = record.atomname;
							if (atmids.count(atm) == 0)
								atmids.insert(std::make_pair(atm, 1));
							else
								atmids[atm]++;
							record.namemodifier = std::to_string(atmids[atm]);
						} // e.g. HOH do not need to change O to O1
						record.chainid = '0';
						record.residueid = record.residueid-1;
						record.atomid = n_aid;
						n_aid++;
						n_records.push_back(record);
					}
				} // every atom
			} // every ligand part
			else // key residues
			{
				for (auto &record : records)
				{
					record.chainid = c + 48;
					record.residueid = r + 1; // s.t. Gly is 0 and n.
				}
				o_records[c-1].push_back(records);
			}
		} // every residue
	} // every chain

	// add gly in the two ends of key residues
	auto rng=NSPdstl::RandomEngine<>::getinstance().realrng(0,1);
	for (auto &oc_rs : o_records)
	{
		auto crd_beg = mccrds(oc_rs[0]);
		auto crd_end = mccrds(oc_rs[oc_rs.size()-1]);
		double phi, psi, b_phi;
		double b_psi = torsion(crd_beg[0], crd_beg[1], crd_beg[3], crd_beg[2])/degree;
		while(true) {
			NSPpdbstatistics::PhiPsiDistr::coildistr().randomphipsi(rng,&phi,&psi);
			while(phi > 180.0)
				phi -= 360.0;
			while(phi < -180.0)
				phi += 360.0;
			while(psi > 180.0)
				psi -= 360.0;
			while(psi < -180.0)
				psi += 360.0;
			if(psi > b_psi - 3.0 && psi < b_psi + 3.0) {
				b_phi = phi;
				break;
			}
		}
		double e_psi;
		if(oc_rs.size() == 1)
			e_psi = b_psi;
		else
			e_psi = torsion(crd_end[0], crd_end[1], crd_end[3], crd_end[2])/degree;
		std::vector<XYZ> crd_pre, crd_bac; // pre_gly, bac_gly
		NSPpdbstatistics::PhiPsiDistr::glydistr().randomphipsi(rng,&phi,&psi);
		genpregly(crd_pre, crd_beg, b_phi, psi);
		NSPpdbstatistics::PhiPsiDistr::glydistr().randomphipsi(rng,&phi,&psi);
		genbackgly(crd_bac, crd_end, e_psi, phi);

		// create record for pre_gly and bac_gly;
		auto preglys = make_pregly(crd_pre, oc_rs[0][0].chainid);
		for (auto &pg : preglys)
		{
			pg.atomid = n_aid;
			n_aid++;
			n_records.push_back(pg);
		}
		for (auto &ocr : oc_rs)
			for (auto &o : ocr)
			{
				o.atomid = n_aid;
				n_aid++;
				n_records.push_back(o);
			}
		auto bacglys = make_bacgly(crd_bac, oc_rs[0][0].chainid, oc_rs.size()+1);
		for (auto &bg : bacglys)
		{
			bg.atomid = n_aid;
			n_aid++;
			n_records.push_back(bg);
		}
	} // every chain (key residue group), add two gly

	// write
	for (auto &r : n_records)
		ofs << r.toString() << std::endl;
}


