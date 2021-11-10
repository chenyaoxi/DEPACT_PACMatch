/*
 * subsite.h
 *
 *  Created on: 2018年8月12日
 *      Author: hyliu
 */

#ifndef SUBSITE_H_
#define SUBSITE_H_
#include "atomcontacts.h"
#include "residuecontactscores.h"
#include "scorecontact.h"
#include "proteinrep/pdbrecord.h"
#include <iostream>
#include <fstream>
using namespace NSPproteinrep;
namespace subsitedesign {
/**
 * @brief A subpocket contains environmental residues that interact with the target ligand
 */
class Subsite {
public:
	/**
	 * @brief a pair-wise interaction between a target atom and an atom from the environment
	 * @note contains atom indices and distance squared
	 */
	struct AtomPairD2 {
		int atom1 { -1 }; /*!< atom idx of the target if smaller than the number of atoms in
		 in targetstruct. Otherwise atom1 is a mediating water or metal ion
		 whose position in keyresidues_ is indicated in mediachain_
		 @see mediachain_
		 */
		int atom2 { -1 }; ///< atom idx in the key residue. always 0 if the key residue is a single atom water or metal ion
		double d2 { 0.0 }; ///< inter-atomic distance squared
		AtomPairD2() {
			;
		}
		AtomPairD2(int a1, int a2, double r2) :
				atom1(a1), atom2(a2), d2(r2) {
			;
		}
	};
	/**
	 * @brief an environmental protein residue or mediating water/ion
	 */
	struct KeyResidue {
		NSPproteinrep::AAConformer conformer; ///< contains atomnames and coordinates
		std::vector<AtomPairD2> apd2s; ///< atom pairwise interactions of this residue with ligand
		std::vector<int> nseparatingbonds; /*!<the smallest number of bonds separating the atom from
		 the atoms contacting the ligand or mediating water.
		 All values above 3 is set to 3.
		 */
	};
	/**
	 * @brief derive subpocket from a template and a substructure aligment
	 * @param targetstruct target
	 * @param acontacts atomic contacts of template with environment
	 * @param al substructure aligment between template and target
	 * @see keyresidues_
	 */
	Subsite(const TargetStruct &targetstruct,
			const std::vector<AtomContacts> & acontacts,
			const TmpltSSAs::Alignment &al);
	/**
	 * @brief merge s1 & s2 and output the lowest merge results. (such merge divided subpocket into pieces).
	 */
	Subsite(const Subsite & s1, const Subsite &s2, std::string merge);

	/**
	 * @brief merge two subpockets into one
	 */
	Subsite(const Subsite & s1, const Subsite &s2);

	// no select, but include all residues as key residues. Used for native pocket
	Subsite(const TargetStruct &targetstruct,
			const std::vector<AtomContacts> & acontacts,
			const TmpltSSAs::Alignment &al, bool all_include);

	/**
	 *@brief returns environmental residues comprising the subpocket
	 *@see keyresidues_
	 */
	std::vector<std::vector<KeyResidue>> &keyresidues() {
		return keyresidues_;
	}
	/**
	 *@brief returns environmental residues comprising the subpocket
	 *@see keyresidues_
	 */
	const std::vector<std::vector<KeyResidue>> &keyresidues() const {
		return keyresidues_;
	}

	ResidueContactScores & keyresiduescores(){
		return keyresiduescores_;
	}

	const ResidueContactScores &keyresiduescores() const {
		return keyresiduescores_;
	}

	double totalscore(double wdmldpl) const
	{
		double sum = 0.0;
		for (auto &pl : keyresiduescores_.plscores)
			sum += keyresiduescores_.totalplscore(pl.first);
		auto &mls = keyresiduescores_.mlscores;
		for (auto iter = mls.begin(); iter != mls.end(); iter++)
		{
			ResidueID resid = iter->first;
			std::string mtype = keyresidues_[resid.first][resid.second].conformer.residuename;
			sum += wdmldpl * keyresiduescores_.totalmlscore(iter->first) *
					ScoreContact::theta_m(mtype, keyresiduescores_.totalmpscore(iter->first));
		}
		return sum;
	}

	std::vector<double> detailtotalscore(std::string outfile, double wdmldpl) const
	{
		std::vector<double> scores; // 0: totalscore; 1: totalplscore; 2: totalpmlscore
		std::ofstream ofs(outfile.c_str(), std::ios::app);
		ofs << "every pl score: " << std::endl;
		double plsum = 0.0;
		for (auto &pl : keyresiduescores_.plscores)
		{
			ofs << pl.first.first << " " << pl.first.second << " " << keyresiduescores_.totalplscore(pl.first) << "; ";
			plsum += keyresiduescores_.totalplscore(pl.first);
		}
		ofs << std::endl;
		ofs << "every pml score:" << std::endl;
		double pmlsum = 0.0;
		auto &mls = keyresiduescores_.mlscores;
		for (auto iter = mls.begin(); iter != mls.end(); iter++)
		{
			ResidueID resid = iter->first;
			std::string mtype = keyresidues_[resid.first][resid.second].conformer.residuename;
			ofs << resid.first << " " << resid.second << " " <<
					wdmldpl * keyresiduescores_.totalmlscore(iter->first) *
					ScoreContact::theta_m(mtype, keyresiduescores_.totalmpscore(iter->first)) << "; ";
			pmlsum += wdmldpl * keyresiduescores_.totalmlscore(iter->first) *
					ScoreContact::theta_m(mtype, keyresiduescores_.totalmpscore(iter->first));
		}
		ofs << std::endl;
		ofs << "plscore: " << plsum << "; pmlscore: " << pmlsum << std::endl;
		scores.push_back(plsum+pmlsum);
		scores.push_back(plsum);
		scores.push_back(pmlsum);
		return scores;
	}

	/**
	 *@brief output the subpocket
	 */
	void writepdb(std::ostream &os) const;
	/**
		 *@brief output the subpocket
		 */
	void writepdb(const std::string &filename) const {
		std::ofstream ofs(filename);
		writepdb(ofs);
	}
	std::vector<AtomContacts> &atomcontacts() const;
	ResidueContactScores &residuecontactscores() const;
	bool goodcontact();

	TargetStruct gettargetstruct() {return *targetstruct_;}

private:
	const TargetStruct* targetstruct_;///! target
	/**
	 * @brief components forming the subpocket
	 * the residues are grouped into pieces (or "chains", as vector of vectors);
		 *each piece in one vector is either a continuous peptide segments of one or more
		 *residues, or a mediating water or ion.
		 *The key residues are selected from the template's atomic contacts
		 *by the constructor.
		 *each protein residue must has certain minimum numbers
		 *of contacts with ligand.
		 *Side chains are removed if there is no sidechain-ligand contact.
		 *Gly linkers are included if two ligand-contacting residues
		 *are separated by no more than two positions in sequence.
		 */
	std::vector<std::vector<KeyResidue>> keyresidues_;
	ResidueContactScores keyresiduescores_;
	///contacts between ligand atoms and atoms in keyresidues_;
//	mutable std::shared_ptr<std::vector<AtomContacts>> atomcontacts_{nullptr};
	///contact scores of key residues
//	mutable std::shared_ptr<ResidueContactScores> residuecontactscores_{nullptr};
	std::vector<int> mediachain_; ///<"chain" idx in keyresidues_ of the mediating water/ion
	//the following functions are used in the constructor

	typedef std::pair<int, int> ResidueID;
	///used in the constructor to select key residues
	bool testresidue(const NSPproteinrep::AAConformersInModel &protein,
			const ResidueID &residueid, const Move3D &move3d,
			std::vector<AtomPairD2> & apd2s, double &score,int &num_np_c_contacts,
			int &num_p_nc_contacts);
	///used in the constructor to select key residues
	std::vector<std::vector<ResidueID>> getpieces(
			const NSPproteinrep::AAConformersInModel &protein,
			const std::set<ResidueID> &aaselected, const Move3D& move3d);
	///stores transformed coordinates of the environmental residues
	std::map<ResidueID, std::vector<NSPgeometry::XYZ>> newcrds_;
	ResidueContactScores ctscores_a_;
	///used in the constructor to select key residues
	void selectresidues(const std::vector<AtomContacts> & acontacts,
			const TmpltSSAs::Alignment &al, std::set<ResidueID>& contactingaa,
			std::map<ResidueID, std::vector<AtomPairD2>> &contactingapd2,
			std::vector<ResidueID> & mediatingresidues);
	///used in the constructor to select key residues: score < 0 is enough to enclude
	void selectresidues_0(const std::vector<AtomContacts> & acontacts,
			const TmpltSSAs::Alignment &al, std::set<ResidueID>& contactingaa,
			std::map<ResidueID, std::vector<AtomPairD2>> &contactingapd2,
			std::vector<ResidueID> & mediatingresidues);
	///used in the constructor to determine bond separation @see KeyResidue
	void setnseparatingbond(std::vector<KeyResidue> &keyresidues);
	double scoremresidue(const NSPproteinrep::AAConformersInModel &protein,
			ResidueID mid,ResidueID prid,const
			std::map<std::pair<int,int>,double> &iatomd2);
};
/**
 * @brief determin if two subsites are compatible
 * Based on whether atoms of different subpockets clash.
 * The clash distances depend on bond separation.
 */
bool compatible(const Subsite &s1, const Subsite &s2);
/**
 * @brief combine subpockets by finding maximal cliques
 * @param ssites are the component subpocket
 * @return vector of cpointers to combined pockets.
 */
std::vector<std::shared_ptr<Subsite>> combinebycliques(
		const std::vector<std::shared_ptr<Subsite>> &ssites);
std::vector<std::shared_ptr<Subsite>> combinebycliques_mm(
		const std::vector<std::shared_ptr<Subsite>> &ssites,
		std::vector<std::vector<int>> &mmbers);
}

#endif /* SUBSITE_H_ */
