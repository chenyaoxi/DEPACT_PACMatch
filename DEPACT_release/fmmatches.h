/*
 * fmmatches.h
 *
 *  Created on: 2018年8月11日
 *      Author: hyliu
 */

#ifndef FMMATCHES_H_
#define FMMATCHES_H_
#include "basicfragment.h"
#include "twowaymatches.h"
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <memory>
namespace myobcode{
/*
 * Find matches between a basic fragment type and a ligand molecule.
 *
 * The basic fragments is a user-supplied SMARTS pattern.
 * @param unique if true and if there are matches containing the same set of atoms,
 * only one is returned (e.g., if the basic fragment is a benzene ring, unique=true
 * will return only one match, unique=false may return 6).
 * When the found matches are to be used to align common substructures in two molecules,
 * Unique should be set to true for one molecule and to false for another molecule,
 * so that permutations of alignment are considered non-redundantly.
 */
typedef TwoWayMatches<subsitedesign::BasicFragment,OpenBabel::OBMol> FMMatches;
std::vector<std::shared_ptr<FMMatches>> findfmmatches(
		const subsitedesign::BasicFragment *sf,OpenBabel::OBMol *mol,bool unique=false);
/*
 * obtain atomtypes for atoms in a molecule based on matches to a set of
 * BasicFragments
 */
//std::vector<std::string> findatomtypes(const std::map<std::string,std::vector<std::shared_ptr<FMMatches>>> &fmmatches);
}




#endif /* FMMATCHES_H_ */
