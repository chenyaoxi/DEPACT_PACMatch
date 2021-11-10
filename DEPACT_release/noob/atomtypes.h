/*
 * atomtypes.h
 *
 *  Created on: 2018年8月16日
 *      Author: hyliu
 */

#ifndef ATOMTYPES_H_
#define ATOMTYPES_H_
#include <vector>
#include <string>
#include <map>
#include <set>
#include <bitset>
#include <iostream>
namespace subsitedesign {
/**
 * @ The detailed atom type,probably defined using SMARTS patterns in an input file
 */
struct AtomType {
	enum Features {
		SPECIFIED,
		ISPOLAR,
		ISNONPOLAR,
		ISHBDONOR,
		ISHBACCEPTOR,
		ISAROMATIC,
		ISNEGATIVE,
		ISPOSITIVE,
		NFEATURES
	};
	typedef std::bitset<NFEATURES> AtomFeatures;
	AtomType(const std::string cnm="",const std::string & nm = "", const std::string &sm="",
			const AtomFeatures &ft = AtomFeatures(0)) :
			codename(cnm),name(nm), smarts(sm),features(ft) {
		;
	}
	/**
	 * @brief returns the atom types read from an input file.
	 *
	 * Employs lazy initialization to read input file
	 */
	static const std::map<std::string,AtomType> & getmap(
			const std::string &filename="atomtypesmarts.txt");
	static void initmap(const std::string &filename,std::map<std::string,AtomType> &map);
	static bool ishydrogen(const std::string &anm) {
		return (anm.substr(0,2)=="H_" || anm.substr(0,2)=="h_");
	}
	std::string name; ///< name of detailed type
	std::string codename0; ///<typename, less specific than name
	std::string codename; ///<name of recoded type. Perhaps more coarsely typed
	std::string smarts; ///<SMARTS string defining the detailed type
	AtomFeatures features { 0 };
	bool featuretest(Features f) const {
		return features.test(f);
	}
	static const AtomType & atomtype(const std::string &nm){
		auto maps = getmap();
		if (maps.count(nm) > 0) return maps.at(nm);
		else
		{
			std::cout << "atomtypesmarts.txt don't have such kind of atom: " << nm << std::endl;
			exit(1);
		}
	}
	void setfeature(Features f) {
		features.set(f);
	}
	void unsetfeature(Features f) {
		features.reset(f);
	}
};
std::string proteinatomtype(const std::string &residue, std::string atom);
std::string proteinatomtype2(const std::string &residue, std::string atom);
std::string ptype2(const std::string & ptype1);
std::string getcodename0(const std::string &atype);
std::string getcodename(const std::string &atype);
}

#endif /* ATOMTYPES_H_ */
