/*
 * subsite.cpp
 *
 *  Created on: 2018年8月12日
 *      Author: hyliu
 */
#include "subsite.h"
#include "basicfragment.h"
#include "atomtypes.h"
#include "scorecontact.h"
#include <algorithm>
#define CLASH2PNC  6.25 ///<  clash dist squared for polar-non-carbon contacts
#define CLASH2PC  8.41 ///< clash dist squared for polar-carbon contacts
#define CLASH2NPNC 8.41 ///< clash dist squared for non-polar-non-carbon contacts
#define CLASH2NPC 9.61 ///< clash dist squared for  non-polar-carbon contacts
#define R2CONTACTP 18.49 ///< contact dist squared for polar target atoms
#define R2CONTACTNP 18.49 //contact dist squared for non-polar target atoms
#define MINNUMPNC 1 ///< minimum number of polar-non-carbon contacts to include a residue
#define MINNUMNPC 3 ///< minimum number of non-polar-carbon contacts to include a residue
#define RESIDUESCORE_CUT -1.0 ///maximum contact score to keep residue
//#define MPSCORE_CUT0 100.1
//#define MPSCORE_CUT1 100.1
using namespace NSPproteinrep;
using namespace subsitedesign;
/**
 * @brief a utility function to obtain transformed coordinates of atoms in an AAConformer
 * @TODO should be placed in another source file for reuse
 */
static std::vector<NSPgeometry::XYZ> newaacrd(const AAConformer &conf,
		const Move3D &move3d) {
	std::vector<NSPgeometry::XYZ> res;
	const std::map<std::string, NSPgeometry::XYZ> & crds = conf.getglobalcrd();
	for (auto & a : conf.atomlist) {
		NSPgeometry::XYZ r0 = crds.at(a);
		res.push_back(move3d.move(r0.tovector()));
	}
	return res;
}
/**
 * @brief determine if atoms are carbon based on atomnames
 */
static std::vector<bool> carbonatoms(
		const std::vector<std::string> &atomnames) {
	std::vector<bool> res;
	for (auto &n : atomnames)
		if (n[0] != 'C')
			res.push_back(false);
		else
			res.push_back(true);
	return res;
}
static void copyresiduescores(const ResidueContactScores &s1,
		ResidueContactScores &s2,
		std::map<std::pair<int,int>,std::pair<int,int>> & idmap){
	for(auto &pl:s1.plscores){
		if(idmap.find(pl.first)==idmap.end()) continue;
		s2.plscores[idmap[pl.first]]=pl.second;
	}
	for(auto &ml:s1.mlscores){
		if(idmap.find(ml.first)==idmap.end()) continue;
		s2.mlscores[idmap[ml.first]]=ml.second;
	}
	for(auto &mp:s1.mpscores){
		if(idmap.find(mp.first)==idmap.end()) continue;
		for(auto &pm:mp.second){
			if(idmap.find(pm.first)==idmap.end()) continue;
			s2.mpscores[idmap[mp.first]][idmap[pm.first]]=pm.second;
		}
	}
}
Subsite::Subsite(const Subsite &s1, const Subsite &s2) {
	assert(s1.targetstruct_ == s2.targetstruct_);
	targetstruct_ = s1.targetstruct_;

	int nkr = s1.keyresidues_.size();  ///< number of pieces forming s1
	for (auto &kvec1 : s1.keyresidues_) {
		keyresidues_.push_back(kvec1);
	}

	mediachain_ = s1.mediachain_;
	for (auto m : s2.mediachain_)
		mediachain_.push_back(m + nkr);
	int ntgatoms = targetstruct_->atomtypes.size(); ///<number of atoms in target
	int nm1 = s1.mediachain_.size();
	for (auto &kvec2 : s2.keyresidues_) {
		keyresidues_.push_back(kvec2);
		auto &kv = keyresidues_.back();
		for (auto &kr : kv) {
			for (auto &apd2 : kr.apd2s) {
				if (apd2.atom1 >= ntgatoms) {
					apd2.atom1 += nm1;  ///<idx of mediating atom from s2 need to be shifted
				}
			}
		}
	}
	keyresiduescores_=s1.keyresiduescores_;
	std::map<ResidueID,ResidueID> idmap;
	for(int c=0;c<s2.keyresidues_.size();++c){
		for(int r=0;r<s2.keyresidues_.at(c).size();++r){
			ResidueID rid(c,r);
			idmap[rid]=ResidueID(c+s1.keyresidues_.size(),r);
		}
	}
	copyresiduescores(s2.keyresiduescores_,keyresiduescores_,idmap);
}

Subsite::Subsite(const TargetStruct &targetstruct,
		const std::vector<AtomContacts> & acontacts,
		const TmpltSSAs::Alignment &al) :
		targetstruct_(&targetstruct) {
	std::set<ResidueID> aaselected; ///<selected residues in template
	std::vector<ResidueID> mediatingresidues; ///<selected mediating water or ion
	std::map<ResidueID, std::vector<AtomPairD2>> contactingapd2;
	selectresidues(acontacts, al, aaselected, contactingapd2,
			mediatingresidues);
	const AAConformersInModel *protein = acontacts[0].model;
	std::vector<std::vector<ResidueID>> pieces = getpieces(*protein, aaselected,
			al.move3d);
//	std::vector<ResidueID> newmediatingids;
	for (auto & mid : mediatingresidues) {
		for (int i = 0; i < pieces.size(); ++i) {
			if (mid == pieces[i][0]) {
				mediachain_.push_back(i);
			}
		}
	}
	//making the key residues
	int nlatom=targetstruct.atomtypes.size();
	std::map<ResidueID,ResidueID> idxmap;

	for (auto&p : pieces) {
		keyresidues_.push_back(std::vector<KeyResidue>());
		auto &krvec = keyresidues_.back();
		for (auto &aa : p) {
			const AAConformer &conf = protein->conformers.at(aa.first).at(
					aa.second);
			if (!conf.isnaturalaa()) {
				if (newcrds_.count(aa) == 0) {
					std::cout << "no such aa in newcrds_" << std::endl;
					continue;
				}//for case like PDB2rcy: aa:{1,11}
			}
			krvec.push_back(KeyResidue());
			KeyResidue &kr = krvec.back();
			ResidueID nid(keyresidues_.size()-1,krvec.size()-1);
			idxmap[aa]=nid;
			if (!conf.isnaturalaa()) {
				kr.conformer = conf;
				kr.conformer.globalcrd[conf.atomlist[0]] = newcrds_.at(aa)[0];
				kr.apd2s = contactingapd2.at(aa);
			} else {
				bool keepsidechain = false;
				auto itr = contactingapd2.find(aa);
				if (itr != contactingapd2.end()) {
					for (auto &apd : itr->second) {
						if(apd.d2>R2CONTACTNP) continue;
						std::string aname = conf.atomlist.at(apd.atom2);
						if (!AAConformer::ismainchain(aname)) {
							keepsidechain = true;
							break;
						}
					}
				}
				if (keepsidechain) {
					kr.conformer = conf;
					std::vector<NSPgeometry::XYZ> &ncrd = newcrds_.at(aa);
					for (int i = 0; i < ncrd.size(); ++i) {
						kr.conformer.globalcrd.at(conf.atomlist[i]) = ncrd[i];
					}
					kr.apd2s = contactingapd2[aa];
				} else {
					kr.conformer.residuename = "GLY";
					kr.conformer.chainid_or = conf.chainid_or;
					kr.conformer.residueid_or = conf.residueid_or;
					kr.conformer.atomlist = AAConformer::mainchainatoms;
					for (auto &anm : AAConformer::mainchainatoms) {
						kr.conformer.globalcrd[anm] = NSPgeometry::XYZ(
								al.move3d.move(
										conf.globalcrd.at(anm).tovector()));
					}
					if (contactingapd2.find(aa) != contactingapd2.end()) {
						kr.apd2s = contactingapd2.at(aa);
					}
				} //keepsidechain
			} //else naturalaa
		} //aa:p
	} //p:pieces
	for (auto &krs : keyresidues_) {
		setnseparatingbond(krs);
	}
	copyresiduescores(ctscores_a_,keyresiduescores_,idxmap);
}

Subsite::Subsite(const TargetStruct &targetstruct,
		const std::vector<AtomContacts> & acontacts,
		const TmpltSSAs::Alignment &al, bool all_include) :
		targetstruct_(&targetstruct) {
	std::set<ResidueID> aaselected; ///<selected residues in template
	std::vector<ResidueID> mediatingresidues; ///<selected mediating water or ion
	std::map<ResidueID, std::vector<AtomPairD2>> contactingapd2;
	if (all_include)
		selectresidues_0(acontacts, al, aaselected, contactingapd2,
			mediatingresidues);
	else
		selectresidues(acontacts, al, aaselected, contactingapd2,
			mediatingresidues);
	const AAConformersInModel *protein = acontacts[0].model;
	std::vector<std::vector<ResidueID>> pieces = getpieces(*protein, aaselected,
			al.move3d);
	std::vector<ResidueID> newmediatingids;
	for (auto & mid : mediatingresidues) {
		for (int i = 0; i < pieces.size(); ++i) {
			if (mid == pieces[i][0]) {
				mediachain_.push_back(i);
			}
		}
	}
	//making the key residues
	int nlatom=targetstruct.atomtypes.size();
	std::map<ResidueID,ResidueID> idxmap;

	for (auto&p : pieces) {
		keyresidues_.push_back(std::vector<KeyResidue>());
		auto &krvec = keyresidues_.back();
		for (auto &aa : p) {
			const AAConformer &conf = protein->conformers.at(aa.first).at(
					aa.second);
			if (!conf.isnaturalaa()) {
				if (newcrds_.count(aa) == 0) {
					std::cout << "no such aa in newcrds_" << std::endl;
					continue;
				}//for case like PDB2rcy: aa:{1,11}
			}
			krvec.push_back(KeyResidue());
			KeyResidue &kr = krvec.back();
			ResidueID nid(keyresidues_.size()-1,krvec.size()-1);
			idxmap[aa]=nid;
			if (!conf.isnaturalaa()) {
				kr.conformer = conf;
				kr.conformer.globalcrd[conf.atomlist[0]] = newcrds_.at(aa)[0];
				kr.apd2s = contactingapd2.at(aa);
			} else {
				bool keepsidechain = false;
				auto itr = contactingapd2.find(aa);
				if (itr != contactingapd2.end()) {
					for (auto &apd : itr->second) {
						if(apd.d2>R2CONTACTNP) continue;
						std::string aname = conf.atomlist.at(apd.atom2);
						if (!AAConformer::ismainchain(aname)) {
							keepsidechain = true;
							break;
						}
					}
				}
				if (keepsidechain) {
					kr.conformer = conf;
					std::vector<NSPgeometry::XYZ> &ncrd = newcrds_.at(aa);
					for (int i = 0; i < ncrd.size(); ++i) {
						kr.conformer.globalcrd.at(conf.atomlist[i]) = ncrd[i];
					}
					kr.apd2s = contactingapd2[aa];
				} else {
					kr.conformer.residuename = "GLY";
					kr.conformer.chainid_or = conf.chainid_or;
					kr.conformer.residueid_or = conf.residueid_or;
					kr.conformer.atomlist = AAConformer::mainchainatoms;
					for (auto &anm : AAConformer::mainchainatoms) {
						kr.conformer.globalcrd[anm] = NSPgeometry::XYZ(
								al.move3d.move(
										conf.globalcrd.at(anm).tovector()));
					}
					if (contactingapd2.find(aa) != contactingapd2.end()) {
						kr.apd2s = contactingapd2.at(aa);
					}
				} //keepsidechain
			} //else naturalaa
		} //aa:p
	} //p:pieces
	for (auto &krs : keyresidues_) {
		setnseparatingbond(krs);
	}
	copyresiduescores(ctscores_a_,keyresiduescores_,idxmap);
}

void Subsite::setnseparatingbond(std::vector<Subsite::KeyResidue> &krs) {
	for (auto &kr : krs) {
		kr.nseparatingbonds.resize(kr.conformer.atomlist.size(), 3);
		for (auto &apd2 : kr.apd2s) {
			if (apd2.atom2 < kr.nseparatingbonds.size())
				kr.nseparatingbonds[apd2.atom2] = 0; //direct contact atoms, ns=0
			//for some cases like PDB1dqa: apd2.atom2 = 4 while kr.nseparatingbonds.size = 4;
		}
	}
	for (auto &kr : krs) {
		for (int i = 0; i < kr.nseparatingbonds.size(); ++i) {
			if (kr.nseparatingbonds[i] == 0)
				continue;
			NSPgeometry::XYZ xi = kr.conformer.globalcrd.at(
					kr.conformer.atomlist.at(i));
			for (auto &kr1 : krs) { //all residues in the same piece
				for (auto &apd2 : kr1.apd2s) {
					if (kr1.conformer.atomlist.size() <= apd2.atom2)
						continue;//for a case PDB4p3q: apd2.atom2 = 4, while size = 4.
					NSPgeometry::XYZ xj = kr1.conformer.globalcrd.at(
							kr1.conformer.atomlist.at(apd2.atom2));
					double r2 = (xi - xj).squarednorm();
					if (r2 <= 3.24)
						kr.nseparatingbonds[i] =
								kr.nseparatingbonds[i] > 1 ?
										1 : kr.nseparatingbonds[i];
					else if (r2 <= 7.29)
						kr.nseparatingbonds[i] =
								kr.nseparatingbonds[i] > 2 ?
										2 : kr.nseparatingbonds[i];
				} //apd2
			} //kr1
		} //i
	} //kr
}
bool Subsite::testresidue(const AAConformersInModel &protein,
		const ResidueID &residueid, const Move3D &move3d,
		std::vector<AtomPairD2> & apd2s, double &score, int &num_np_c, int &num_p_nc) {
	const AAConformer &target = targetstruct_->conformer; ///< contains target atom names and coordinates
	const std::vector<std::string> &tatoms = target.atomlist;
	std::vector<bool> tpolar; ///<target atom polar?
	for (auto &n : targetstruct_->atomtypes)
		tpolar.push_back(
				AtomType::atomtype(n).featuretest(
						AtomType::ISPOLAR));
	const AAConformer &residue = protein.conformers.at(residueid.first).at(
			residueid.second);
	if (!residue.crdcomplete()) ///<ignore residues with incomplete coordinates
		return false;
	if (newcrds_.find(residueid) == newcrds_.end()) { ///<lazy update of transformed coordinates
		newcrds_[residueid] = newaacrd(residue, move3d);
	}
	std::vector<NSPgeometry::XYZ> & newcrd = newcrds_[residueid];
	std::vector<bool> carbons = carbonatoms(residue.atomlist); ///<template atom is carbon?

	bool clashed = false;
	score=0.0;
	num_np_c = 0;
	num_p_nc = 0;
	std::string rname=residue.residuename;
	///check clash with target
	for (int i = 0; i < tatoms.size(); ++i) {
		NSPgeometry::XYZ xi = target.getglobalcrd().at(tatoms.at(i));
		std::string tatype=targetstruct_->atomtypes.at(i);
		std::string paname=rname;
		std::string patype=rname;
		for (int j = 0; j < newcrd.size(); ++j) {
			NSPgeometry::XYZ xj = newcrd[j];
			double rcut2 = CLASH2NPNC;
			if (tpolar[i]) {
				if (carbons[j])
					rcut2 = CLASH2PC;
				else
					rcut2 = CLASH2PNC;
			} else if (carbons[j])
				rcut2 = CLASH2NPC;
			double r2 = (xi - xj).squarednorm();
			if (r2 < rcut2)
			{
				if (residue.isnaturalaa()) // as metal coordination is quite close
				{
					clashed = true;
					return false;
				}
				else // metal shouldn't be too close to Carbon
				{
					if (tatoms[i][0] == 'C')
					{
						clashed = true;
						return false;
					}
				}
			}

			// change..
			rcut2 = 25; //tpolar[i] ? R2CONTACTP : R2CONTACTNP;
			if (r2 < rcut2) {
				apd2s.push_back(AtomPairD2(i, j, r2));
				if (tpolar[i] && !carbons[j])
					num_p_nc++;
				if (!tpolar[i] && carbons[j])
					num_np_c++;
 			}
			double r=sqrt(r2);
			double s;
			if(residue.isnaturalaa()){
				paname=residue.atomlist.at(j);
				patype=proteinatomtype(rname,paname);
				s=ScoreContact::score(tatype,patype,r);
				ctscores_a_.insertplscore(residueid,paname,tpolar.size(),i,s);
			} else {
				s=ScoreContact::score_scaled(tatype,patype,r);
				ctscores_a_.insertmlscore(residueid,tpolar.size(),i,s);
			}
			score +=s;
		} //newcrd
	} //tatoms

	return true;
}

std::vector<std::vector<Subsite::ResidueID>> Subsite::getpieces(
		const NSPproteinrep::AAConformersInModel &protein,
		const std::set<ResidueID> &aaselected, const Move3D& move3d) {
	std::vector<std::vector<ResidueID>> pieces;
	///sort the selected residues, to find consecutive pieces
	std::vector<ResidueID> sortedselected;
	for (auto &aa : aaselected)
		sortedselected.push_back(aa);
	std::sort(sortedselected.begin(), sortedselected.end());
	bool newpiece = false;
	ResidueID prevr(-1, -1);
	const AAConformer &target = targetstruct_->conformer;
	const std::vector<std::string> &tatoms = target.atomlist;
	for (int i = 0; i < sortedselected.size(); ++i) {
		ResidueID & rid = sortedselected[i];
		const AAConformer &conformer = protein.conformers[rid.first][rid.second];
		if (i == 0 || !conformer.isnaturalaa() || rid.first != prevr.first
				|| (rid.second - prevr.second) > 3) {
			newpiece = true;
		}
		if (!newpiece) {
			const AAConformer *cn = &conformer;
			//check peptide connection
			for (int r = rid.second - 1; r >= prevr.second; --r) {
				const AAConformer &cfr = protein.conformers[rid.first][r];
				if (!cfr.mainchaincrdcomplete() || !cfr.connectedto(*cn)) {
					newpiece = true;
					break;
				}
				cn = &cfr;
			}
			if (!newpiece) {
				//check steric clash between connecting backbone and the target
				for (int r = prevr.second + 1; r < rid.second; ++r) {
					std::vector<NSPgeometry::XYZ> bcrd = protein.conformers.at(
							rid.first).at(r).getbackbonecrd();
					for (auto &x : bcrd) {
						x = NSPgeometry::XYZ(move3d.move(x.tovector()));
					}
					for (int i = 0; i < tatoms.size(); ++i) {
						NSPgeometry::XYZ xi = target.getglobalcrd().at(
								tatoms.at(i));
						for (int j = 0; j < bcrd.size(); ++j) {
							NSPgeometry::XYZ xj = bcrd[j];
							double rcut2 = CLASH2NPC;
							double r2 = (xi - xj).squarednorm();
							if (r2 < rcut2) {
								newpiece = true;
								break;
							}
						}
						if (newpiece)
							break;
					}
					if (newpiece)
						break;
				}
			}
		} // done deciding on if to start new piece
		if (newpiece) {
			pieces.push_back(std::vector<ResidueID>());
			pieces.back().push_back(rid);
			newpiece = false;
		} else {
			auto &piece = pieces.back();
			for (int r = prevr.second + 1; r <= rid.second; ++r) {
				piece.push_back(std::make_pair(rid.first, r));
			}
		}
		prevr = rid;
	}
	return pieces;
}
void Subsite::selectresidues(const std::vector<AtomContacts> & acontacts,
		const TmpltSSAs::Alignment &al, std::set<ResidueID>& contactingaa,
		std::map<ResidueID, std::vector<AtomPairD2>> &contactingapd2,
		std::vector<ResidueID> &mediatingresidue) {
	AtomContacts merged = mergeatomcontacts(acontacts, al);
	const AAConformersInModel *protein = acontacts[0].model;
	for (auto &c : merged.dcontacts) {
		int ic = c.ichain;
		int ir = c.iresidue;
		ResidueID residueid(ic, ir);
		std::vector<AtomPairD2> apd2s;
		int num_p_nc;
		int num_np_c;
		double score;
		if (testresidue(*protein, residueid, al.move3d, apd2s, score, num_np_c,
				num_p_nc)) {
			if (score<RESIDUESCORE_CUT) //&& (num_np_c >= MINNUMNPC || num_p_nc >= MINNUMPNC))
			{
//				std::cout << score << " " << ic << " " << std::to_string(ir) << std::endl;
				contactingaa.insert(residueid);
				contactingapd2[residueid] = apd2s;
			}
		}
	}
	int widx = targetstruct_->atomtypes.size();

	// for water-media not be erased
	std::set<ResidueID> mmids;

	for (auto &cm : merged.mcontacts) {
		if (contactingaa.find(cm.first) == contactingaa.end())
			continue;
		int mc = cm.first.first;
		int mr = cm.first.second;
		ResidueID mid(mc, mr);
		std::string mtype=protein->conformers.at(mc).at(mr).residuename;
		//transformed coordinate of mediating atom
		NSPgeometry::XYZ xm = newcrds_.at(mid)[0];
		bool keep = true;
		std::map<ResidueID, std::vector<AtomPairD2>> mcapd2s;
		double mpscore=0.0;
		for (auto &c : cm.second) { //loop over indirect contacting residues
			int ic = c.ichain;
			int ir = c.iresidue;
			ResidueID residueid(ic, ir);
//			if (contactingaa.find(residueid) == contactingaa.end()) {
			{
				std::vector<AtomPairD2> apd2s;
				int num_p_nc;
				int num_np_c;
				double score;
				if (testresidue(*protein, residueid, al.move3d, apd2s,
						score,num_np_c,
						num_p_nc)) {
//					contactingaa.insert(residueid);
//					contactingapd2[residueid] = apd2s;
					double s=scoremresidue(*protein,mid,residueid,c.iatomd2);
//					std::cout << "scoremresidue " << s << std::endl;
					if (s == 0) // for water-metal complex
					{
						const std::vector<std::vector<AAConformer>> &aas = protein->conformers;
						const NSPproteinrep::AAConformer &aa=
									aas.at(residueid.first).at(residueid.second);
						std::string rname = aa.residuename;
						if (rname == "HOH")
						{
							mmids.insert(residueid);
							mcapd2s[residueid] = apd2s;
						}
					}
					else if(contactingaa.find(residueid) != contactingaa.end())
					{
						mcapd2s[residueid] = apd2s;
						mpscore+=s;
					}
					else if(s<=ScoreContact::mp_cutoffs(mtype))
					{
						double r_min = 10.0;
						for (auto &apd2 : apd2s)
							if (sqrt(apd2.d2) < r_min)
								r_min = sqrt(apd2.d2);
						if (r_min > 3.5) continue;
						mcapd2s[residueid] = apd2s;
						mpscore+=s;
					}
				}
			}
		}
//		}
		if(ScoreContact::theta_m(mtype,mpscore)<0.5) keep=false; // ori. is 0.5
		if (keep) {
			for (auto &c : cm.second) { //loop over indirect contacting residues
				int ic = c.ichain;
				int ir = c.iresidue;
				ResidueID residueid(ic, ir);
				if(mcapd2s.find(residueid)==mcapd2s.end()) continue;
				if (contactingaa.find(residueid) == contactingaa.end()) {
					contactingaa.insert(residueid);
					contactingapd2[residueid] = mcapd2s[residueid];
				}
				for (auto & ijd2 : c.iatomd2) {
					assert(ijd2.first.first == 0);
					contactingapd2[residueid].push_back(
							AtomPairD2(widx, ijd2.first.second, ijd2.second));
				}
			}
			mediatingresidue.push_back(mid);
			widx++;
		} else {
			if (mmids.count(mid) != 0) continue;
			contactingaa.erase(mid);
			contactingapd2.erase(mid);
		} //keep
	}	 //mcontacts
}

void Subsite::selectresidues_0(const std::vector<AtomContacts> & acontacts,
		const TmpltSSAs::Alignment &al, std::set<ResidueID>& contactingaa,
		std::map<ResidueID, std::vector<AtomPairD2>> &contactingapd2,
		std::vector<ResidueID> &mediatingresidue) {
	AtomContacts merged = mergeatomcontacts(acontacts, al);
	const AAConformersInModel *protein = acontacts[0].model;
	for (auto &c : merged.dcontacts) {
		int ic = c.ichain;
		int ir = c.iresidue;
		ResidueID residueid(ic, ir);
		std::vector<AtomPairD2> apd2s;
		int num_p_nc;
		int num_np_c;
		double score;
		if (testresidue(*protein, residueid, al.move3d, apd2s, score, num_np_c,
				num_p_nc)) {
			// changed...
//			if (score < 0.0)
//			{
				contactingaa.insert(residueid);
				contactingapd2[residueid] = apd2s;
//			}
		}
	}
	int widx = targetstruct_->atomtypes.size();

	// for water-media not be erased
	std::set<ResidueID> mmids;

	for (auto &cm : merged.mcontacts) {
		if (contactingaa.find(cm.first) == contactingaa.end())
			continue;
		int mc = cm.first.first;
		int mr = cm.first.second;
		ResidueID mid(mc, mr);
		std::string mtype=protein->conformers.at(mc).at(mr).residuename;
		//transformed coordinate of mediating atom
		NSPgeometry::XYZ xm = newcrds_.at(mid)[0];
		bool keep = true;
		std::map<ResidueID, std::vector<AtomPairD2>> mcapd2s;
		double mpscore=0.0;
		for (auto &c : cm.second) { //loop over indirect contacting residues
			int ic = c.ichain;
			int ir = c.iresidue;
			ResidueID residueid(ic, ir);
//			if (contactingaa.find(residueid) == contactingaa.end()) {
			{
				std::vector<AtomPairD2> apd2s;
				int num_p_nc;
				int num_np_c;
				double score;
				if (testresidue(*protein, residueid, al.move3d, apd2s,
						score,num_np_c,
						num_p_nc)) {
//					contactingaa.insert(residueid);
//					contactingapd2[residueid] = apd2s;
					double s=scoremresidue(*protein,mid,residueid,c.iatomd2);
//					std::cout << "scoremresidue " << s << std::endl;
					if (s == 0) // for water-metal complex
					{
						const std::vector<std::vector<AAConformer>> &aas = protein->conformers;
						const NSPproteinrep::AAConformer &aa=
									aas.at(residueid.first).at(residueid.second);
						std::string rname = aa.residuename;
						if (rname == "HOH")
						{
							mmids.insert(residueid);
							mcapd2s[residueid] = apd2s;
						}
					}
//					std::cout << "scoremresidue " << s << std::endl;
					else if(contactingaa.find(residueid) != contactingaa.end())
					{
						mcapd2s[residueid] = apd2s;
						mpscore+=s;
					}
//					else if(s<=ScoreContact::mp_cutoffs(mtype))
					else
					{
						double r_min = 10.0;
						for (auto &apd2 : apd2s)
							if (sqrt(apd2.d2) < r_min)
								r_min = sqrt(apd2.d2);
						if (r_min > 3.5) continue;
						mcapd2s[residueid] = apd2s;
						mpscore+=s;
					}
				}
			}
		}
//		if(mpscore > 0.0) keep=false;
		if (keep) {
			for (auto &c : cm.second) { //loop over indirect contacting residues
				int ic = c.ichain;
				int ir = c.iresidue;
				ResidueID residueid(ic, ir);
				if(mcapd2s.find(residueid)==mcapd2s.end()) continue;
				if (contactingaa.find(residueid) == contactingaa.end()) {
					contactingaa.insert(residueid);
					contactingapd2[residueid] = mcapd2s[residueid];
				}
				for (auto & ijd2 : c.iatomd2) {
					assert(ijd2.first.first == 0);
					contactingapd2[residueid].push_back(
							AtomPairD2(widx, ijd2.first.second, ijd2.second));
				}
			}
			mediatingresidue.push_back(mid);
			widx++;
		} else {
			if (mmids.count(mid) != 0) continue;
			contactingaa.erase(mid);
			contactingapd2.erase(mid);
		} //keep
	}	 //mcontacts
}
void Subsite::writepdb(std::ostream &os) const {
	const AAConformer &target = targetstruct_->conformer;
	int resid = 1;
	int atomid = 1;
	char insertioncode = ' ';
	std::vector<char> chainids {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
		'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
		'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};
	int cid = 0;
	std::vector<PdbRecord> rcrdstg = target.make_pdbrecords(chainids[cid],
			resid, insertioncode, atomid);
	for (auto &rc : rcrdstg) {
		os << rc.toString() << std::endl;
	}
	atomid += target.atomlist.size();
	for (auto &kv : keyresidues_) {
		cid++;
		resid = 1;
		for (auto &kr : kv) {
			std::vector<PdbRecord> rcrds = kr.conformer.make_pdbrecords(
					chainids[cid], resid++, insertioncode, atomid);
			atomid += kr.conformer.atomlist.size();
			for (auto &rc : rcrds) {
				os << rc.toString() << std::endl;
			}
		}
	}
}
double Subsite::scoremresidue(const NSPproteinrep::AAConformersInModel &protein,
		ResidueID mid,ResidueID prid,
		const std::map<std::pair<int,int>,double> &iatomd2){
	const std::vector<std::vector<AAConformer>> &aas=protein.conformers;
	std::string matype=aas.at(mid.first).at(mid.second).residuename;
	const NSPproteinrep::AAConformer &aa=
				aas.at(prid.first).at(prid.second);
	std::string rname=aa.residuename;
	double score=0.0;
	for(auto &ad2:iatomd2){
			std::string paname=aa.
						atomlist.at(ad2.first.second);
			double r=sqrt(ad2.second);
			if (rname == "HOH" && r < 2.5) // water-metal complex in 2.5 angstrom
			{
				ctscores_a_.insertmpscore(mid, prid, paname, 0.0);
				continue;
			}
			std::string patype=proteinatomtype(rname,paname);
			double s=ScoreContact::score(matype,patype,r);
			ctscores_a_.insertmpscore(mid,prid,paname,s);
			score +=s;
		}
	return score;
}

bool Subsite::goodcontact()
{
	bool result = true;
	auto &tar_conf = targetstruct_->conformer;
	auto &tar_crds = tar_conf.getglobalcrd();
	for (auto &krs : keyresidues_) // chain
	{
		auto &krb_conf = krs[0].conformer;
		auto &kre_conf = krs[krs.size()-1].conformer;
		if (krb_conf.atomlist.size() == 1 || kre_conf.atomlist.size() == 1)
			continue;
		auto &krb_crds = krb_conf.getglobalcrd();
		auto &kre_crds = kre_conf.getglobalcrd();
		double d_b = 4.5;
		double d_e = 4.5;
		std::vector<std::string> min_b(2), min_e(2);
		for (auto iter_t = tar_crds.begin(); iter_t != tar_crds.end(); iter_t++)
		{
			for (auto iter_b = krb_crds.begin(); iter_b != krb_crds.end(); iter_b++)
			{
				double d = iter_t->second.distance(iter_b->second);
				if (d < d_b)
				{
					d_b = d;
					min_b[0] = iter_t->first;
					min_b[1] = iter_b->first;
				}
			}
			for (auto iter_e = kre_crds.begin(); iter_e != kre_crds.end(); iter_e++)
			{
				double d = iter_t->second.distance(iter_e->second);
				if (d < d_e)
				{
					d_e = d;
					min_e[0] = iter_t->first;
					min_e[1] = iter_e->first;
				}
			}
		}
		if (d_b == 4.5)
		{
			std::cout << "bad contact in residue pos: " << krb_conf.chainid_or << " "
					<< krb_conf.residueid_or << std::endl;
			result = false;
		}
		if (d_e == 4.5)
		{
			std::cout << "bad contact in residue pos: " << kre_conf.chainid_or << " "
					<< kre_conf.residueid_or << std::endl;
			result = false;
		}
	}
	return result;
}

Subsite::Subsite(const Subsite &s1, const Subsite &s2, std::string merge)
{
	assert(s1.targetstruct_ == s2.targetstruct_);
	targetstruct_ = s1.targetstruct_;
	int ntgatoms = targetstruct_->atomtypes.size(); ///<number of atoms in target
	int nkr = s1.keyresidues_.size();  ///< number of pieces forming s1
	for (auto &kvec1 : s1.keyresidues_) {
		keyresidues_.push_back(kvec1);
	}
	mediachain_ = s1.mediachain_;
	int m_id = mediachain_.size();
	std::set<int> discard_s2_chain;
	// mediachain has clash?
	for (auto m2 : s2.mediachain_)
	{
		//bool clash = false;
//		bool clash = true;
//		for (auto m1 : s1.mediachain_)
//		{
//			auto kvec1 = s1.keyresidues_[m1];
//			auto kvec2 = s2.keyresidues_[m2];
//			for (auto k1 : kvec1)
//				for (auto k2 : kvec2)
//				{
//					auto c1 = k1.conformer;
//					auto c2 = k2.conformer;
//					for (auto i1 = c1.globalcrd.begin(); i1 != c1.globalcrd.end(); i1++)
//						for (auto i2 = c2.globalcrd.begin(); i2 != c2.globalcrd.end(); i2++)
//							if (i1->second.distance(i2->second) < 3.5)
//								clash = true;
//					if (clash) break;
//				}
//			if (clash) break;
//		}
//		if (!clash)
//		{
//			mediachain_.push_back(nkr++);
//			keyresidues_.push_back(s2.keyresidues_[m2]);
//			auto &kv = keyresidues_.back();
//			for (auto &kr : kv)
//				for (auto &apd2 : kr.apd2s)
//					if (apd2.atom1 >= ntgatoms)
//						apd2.atom1 += m_id;
//		}
//		else discard_s2_chain.insert(m2);
		discard_s2_chain.insert(m2);
	}
	// direct contacting residues.
	for (int j = 0; j < s2.keyresidues_.size(); j++)
	{
		bool med = false;
		for (auto m : s2.mediachain_)
			if (j == m) med = true;
		if (!med)
		{
			auto &kvec2 = s2.keyresidues_[j];
			bool clash = false;
			for (auto kvec1 : keyresidues_)
			{
				for (auto kr1 : kvec1)
					for (auto kr2 : kvec2)
					{
						for (auto i1 = kr1.conformer.globalcrd.begin(); i1 != kr1.conformer.globalcrd.end(); i1++)
							for (auto i2 = kr2.conformer.globalcrd.begin(); i2 != kr2.conformer.globalcrd.end(); i2++)
								if (i1->second.distance(i2->second) < 3.5)
									clash = true;
					}
				if (clash) break;
			}
			if (!clash)
				keyresidues_.push_back(kvec2);
			else discard_s2_chain.insert(j);
		}
	}
	// score copy.
	keyresiduescores_=s1.keyresiduescores_;
	std::map<ResidueID,ResidueID> idmap;
	for(int c=0;c<s2.keyresidues_.size();++c){
		if (discard_s2_chain.count(c) > 0) continue;
		for(int r=0;r<s2.keyresidues_.at(c).size();++r){
			ResidueID rid(c,r);
			idmap[rid]=ResidueID(c+s1.keyresidues_.size(),r);
		}
	}
//	copyresiduescores(s2.keyresiduescores_,keyresiduescores_,idmap);
	for(auto &pl:s2.keyresiduescores_.plscores){
		if(idmap.find(pl.first)==idmap.end()) continue;
		keyresiduescores_.plscores[idmap[pl.first]]=pl.second;
	}

}
