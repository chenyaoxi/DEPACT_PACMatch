/*
 * RotamerLib.h
 *
 *  Created on: 2017��10��23��
 *      Author: notxp
 */

#ifndef DESIGNSEQ_ROTAMERLIB_H_
#define DESIGNSEQ_ROTAMERLIB_H_

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include "designseq/Rotamer.h"
#include "designseq/RotamerGroup.h"
#include "designseq/RotamerEnergyBBDep.h"
#include "designseq/ResName.h"
#include "designseq/ProteinRep.h"
#include "dataio/datapaths.h"


namespace NSPdesignseq {

using namespace std;

class RotamerLib {
public:
	vector<RotamerGroup*> aaRotGroups;

	RotamerGroup allRots;
	map<string,RotamerEnergyBBDep*> rotEnergyMap;
	map<string, int> rotNameToIndex;
	ResName rn;

	RotamerLib();
	RotamerLib(const string& rotLibType);
	/*
	 * Read and store the required rotamer lib on first call.
	 * returns the read lib on later calls.
	 */
   static RotamerLib & rotamerlib(const string &rotLibType);
   /*
    * choose a rotamer lib based on the backbone torsional angles
    */
	static RotamerLib &rotamerlibpp(double phi,
			double psi){
		Phipsi pp(phi,psi);
		std::string rotlibtype="X1";
		rotlibtype[0] = pp.regionAB();
		return rotamerlib(rotlibtype);
	}
	RotamerGroup* getAAGroup(string& triName);
	float getRotamerEnergy(string& rotName, int ppType) const;

	virtual ~RotamerLib();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_ROTAMERLIB_H_ */
