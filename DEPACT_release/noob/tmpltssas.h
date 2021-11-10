/*
 * tmpltssas.h
 *
 *  Created on: 2018年8月9日
 *      Author: hyliu
 */

#ifndef TMPLTSSAS_H_
#define TMPLTSSAS_H_
#include "move3d.h"
#include "proteinrep/aaconformer.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <string>
#include <map>
#include <set>
#include <vector>

namespace subsitedesign{
/**
 * @brief info about a target ligand
 *
 */
struct TargetStruct{
	TargetStruct(){;}
	TargetStruct(std::string name,
			const std::vector<std::string> & atypes,std::istream &is);
	std::string molname; ///< an identifier
	std::vector<std::string> atomtypes; ///<name of atom types
	NSPproteinrep::AAConformer conformer; ///< contains atomnames and coordinates
	template<class Archive>
	    void serialize(Archive & ar, const unsigned int version)
	    {
	        ar & molname;
	        ar & atomtypes;
	        ar & conformer;
	    }
};
/*
 * Substructure alignments in a template molecule to a given target molecule
 */
struct TmpltSSAs{
	std::string tmpltname; //molecule title of the template, follow PDB chemical component dictionary convention
	std::string targetname;
	std::vector<std::string> atomtypes;//atom types in the template moleucule, as defined in basicfragments.txt
	struct Alignment{
		std::set<std::pair<int,int>> alignedpairs;
		Move3D move3d;
		double rmsd;
		template<class Archive>
		    void serialize(Archive & ar, const unsigned int version)
		    {
		        ar & alignedpairs;
		        ar & move3d;
		        ar & rmsd;
		    }
	};
	std::vector<Alignment> alignments;
	enum namefields {PDBID,COMPONENTID,MODELNO,CHAINID,RESIDUENO};
	static std::vector<std::string> parsetmpltname(const std::string &name);
	template<class Archive>
	    void serialize(Archive & ar, const unsigned int version)
	    {
	        ar & tmpltname;
	        ar & targetname;
	        ar & atomtypes;
	        ar & alignments;
	    }
};
}



#endif /* TMPLTSSAS_H_ */
