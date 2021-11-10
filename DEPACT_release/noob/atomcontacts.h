/*
 * subsite.h
 *
 *  Created on: 2018年8月10日
 *      Author: hyliu
 */

#ifndef ATOMCONTACTS_H_
#define ATOMCONTACTS_H_
#include "proteinrep/aaconformer.h"
#include "tmpltssas.h"
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>

namespace subsitedesign{
/**
 * @brief describes contacts between one or more atoms in a ligand
 * with environment in a protein-liand complex.
 *
 * Indirect contacts (e.g, mediated through a water or a metal ion) is included.
 * Any non-protein residue contains exactly one atom is considered a potential
 * mediating atom
 */
struct AtomContacts {
	/**
	 * @brief describes contacts between one or more atoms in a ligand
	 * with possibly multiple atoms in one environmental residue
	 */
	struct Contact{
		int ichain{-1};  ///< chain number of the pdb residue
		int iresidue{-1}; ///<residue number of the pdb residue
		int lchain{-1}; ///<chain number of the focused ligand residue
		int lresidue{-1}; ///<residue number of the focused ligand residue
		std::map<std::pair<int,int>,double> iatomd2; /*!< atom indices in the
						contacting ligand and protein residue and distance squared*/
		Contact(){;}
		Contact(int lc,int lr):lchain(lc),lresidue(lr){;}
		bool operator<(const Contact &c) const { ///<for sorting based on protein chain and residue numbers
			if(ichain<c.ichain) return true;
			else if(ichain==c.ichain){
				if(iresidue <c.iresidue) return true;
			}
			return false;
		}
		bool sameresidue(const Contact &c) const { ///<if c refers to the same environmental residue as this
			return(ichain==c.ichain && iresidue==c.iresidue);
		}
		/**
		 * @brief for archiving using boost library
		 */
		template<class Archive>
		    void serialize(Archive & ar, const unsigned int version)
		    {
		        ar & ichain;
		        ar & iresidue;
		        ar & lchain;
		        ar & lresidue;
		        ar & iatomd2;
		    }
	};
	const NSPproteinrep::AAConformersInModel *model;  ///<the complete PDB model, usually read from a PDB file
	std::vector<Contact> dcontacts; ///<direct contacts
	std::map<std::pair<int,int>,std::vector<Contact>> mcontacts; /*!<contacts mediated by an one-atom entity, e.g., water
	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 The key is the chain and residue id of the mediating entity*/
	/**
	 * @brief for archiving using boost library
	 */
	template<class Archive>
	void serialize(Archive &ar,const unsigned int version){
		ar &dcontacts;
		ar &mcontacts;
	}
	/**
	 * @brief number of direct contacts
	 */
	int ndcontacts() const {return dcontacts.size();}
	/**
	 * @brief number of indirect contacts
	 */
	int nmcontacts() const {return mcontacts.size();}
	AtomContacts(){;}
	/**
	 * @brief main constructor
	 * @param mdl the PDB model
	 * @param focuschain chain number of the ligand in PDB
	 * @param focusresidue residue number of ligand in PDB
	 * @param atom index of the ligand atom to consider
	 */
	AtomContacts(const NSPproteinrep::AAConformersInModel *mdl,int focuschain,int focusresidue,
			int focusatom);
};
/**
 * @details of a protein atom-ligand atom contact
 *
 * used for storing pre-calculated contacts and for carrying out statistical analysis
 */
struct ContactDetails {
	std::string pdbid;
	int lcid, lrid, laid; ///<ligand chainid,residue id,atom id;
	std::string latype; ///< ligand atom typename
	int indirect; ///< 0 for direct, 1 for indirect;
	int pcid, prid; ///<protein residue ids;
	std::string prname; ///<protein residue name
	std::string paname; ///<protein atom name
	double distance; ///<interatomic distance, protein atom -(mediating) ligand atom
	//the following members only meaningful for mediated contact
	int mcid, mrid; ///<mediating ligand cid,rid;
	std::string mligand; ///<mediating ligand name;
	double mdistance; ///<interatomic distance, ligand-mediating atom
	std::string tostring() const{
		std::string res;
		res=pdbid+" "+std::to_string(lcid)+" "+std::to_string(lrid)+
				" "+std::to_string(laid)+" "+latype+" "+std::to_string(indirect)+
				" " +std::to_string(pcid)+" "+std::to_string(prid)+" "+
				prname +" "+paname+ " "+ std::to_string(distance);
		if(indirect==1){
			res=res+" "+std::to_string(mcid)+" "+std::to_string(mrid)+
					" "+mligand+" "+std::to_string(mdistance);
		}
		return res;
	}
	void read(std::istream &is){
		is>>pdbid >>lcid >>lrid>>laid>>latype>>indirect>>pcid>>prid>>prname
		 >>paname>>distance;
		if(indirect==1)
			is>>mcid>>mrid>>mligand>>mdistance;
	}
};
/**
 * @brief collective contains AtomContacts for all atoms in a template ligand
 */
struct TmpltContacts{
	std::string tmpltname; ///<chemical-component-dictionary name of the template
	std::shared_ptr<std::vector<AtomContacts>> atomcontacts; ///< AtomContacts of all template atoms
	std::shared_ptr<NSPproteinrep::AAConformersInModel> model; ///<PDB model, read from a PDB file
	void readpdb(bool calccontacts=true);
	TmpltContacts(){;}
	TmpltContacts(const std::string &tnm):tmpltname(tnm){
		readpdb(true);
	}
	bool good() const {
		return model != nullptr;
	}
	///for archive with boost library
	BOOST_SERIALIZATION_SPLIT_MEMBER();
	template<class Archive>
	 void save(Archive &ar,const unsigned int version)const {
		ar &tmpltname;
		ar & *atomcontacts;
	}
	template<class Archive>
	 void load(Archive &ar,const unsigned int version){
		ar &tmpltname;
		readpdb(false);
		if(!model) {
			atomcontacts=nullptr;
			return;
		}
		atomcontacts=std::shared_ptr<std::vector<AtomContacts>>
				(new std::vector<AtomContacts>());
		ar & *atomcontacts;
		for(auto &a:*atomcontacts){
			a.model=model.get();
		}
	}
};
/*
 * @brief find AtomContacts for all atoms.
 * mdl is the protein structure
 * ichain and iresidue is the chain number of residue id of the focus residue
 *
 */
std::shared_ptr<std::vector<AtomContacts>> findcontacts(
		const NSPproteinrep::AAConformersInModel *mdl,int ichain,int iresidue);
std::shared_ptr<std::vector<AtomContacts>> findcontacts(
		const NSPproteinrep::AAConformersInModel *mdl,int ichain,int iresidue,
		std::vector<int> selectedatms);
/**
 * @brief merge the AtomContacts of different ligand atoms into one AtomContacts object
 *
 * The same environmental residues contacting different ligand atoms are merged
 *
 */
AtomContacts mergeatomcontacts(const std::vector<AtomContacts> & acontacts,
		const std::vector<int> &atoms);

/**
 * merge AtomContacts for atoms included in a substructure alignment.
 * The original ligand atom information in Contact will be preserved.
 */
inline AtomContacts mergeatomcontacts(const std::vector<AtomContacts> & acontacts,
		const TmpltSSAs::Alignment &al){
		std::vector<int> alnatoms;
		for(auto &aa:al.alignedpairs) alnatoms.push_back(aa.second);
		return mergeatomcontacts(acontacts,alnatoms);
}

/**
 * find ligand atoms colvalently bonded to protein
 */
std::set<int> bondedligandatoms(const std::vector<ContactDetails> &dtls);
/**
 * find atoms covalently connected to the protein and their first and second bonded neighbors
 * @param con is the connectivity of the ligand molecule
 */
std::set<int> extbondedligandatoms(const std::vector<ContactDetails> &dtls,
		const std::vector<std::vector<int>> &con);
std::vector<ContactDetails> collectdetails(const TmpltContacts &tc,
		const std::vector<std::vector<std::string>> &atypes);
double scorelpcontact(const ContactDetails &dtl);
double scorelmcontact(const ContactDetails &dtl);
std::vector<ContactDetails> getdtlsnextmol(std::ifstream &ifs, std::string &molname, int &nlatoms,
		std::vector<std::string> &atypes,std::set<int> &excludedatoms);
}



#endif /* ATOMCONTACTS_H_ */
