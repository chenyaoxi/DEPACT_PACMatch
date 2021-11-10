/*
 * atomtypessmarts.h
 *
 *  Created on: 2018年8月16日
 *      Author: hyliu
 */

#ifndef ATOMTYPESSMARTS_H_
#define ATOMTYPESSMARTS_H_
#include "atomtypes.h"
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
namespace myobcode{
/**
 * @brief To determine the atom types in a ligand based on user-supplied SMARTS patterns defined
 * in "atomtypesmarts.txt"
 * @return atom types of every atom in the molecule.
 * If results[i] contains one element,a unique atom type for atom i can be determined.
 * If results[i] is empty, the atom type of atom i cannot be determined.
 */
std::vector<std::vector<std::string>> findatomtypes(OpenBabel::OBMol *mol);
std::vector<std::vector<int>> connectivity(OpenBabel::OBMol *mol);
}


#endif /* ATOMTYPESSMARTS_H_ */
