/*
 * RotamerEnergyBBDep.cpp
 *
 *  Created on: 2017��10��23��
 *      Author: notxp
 */

#include "designseq/RotamerEnergyBBDep.h"

namespace NSPdesignseq {

RotamerEnergyBBDep::RotamerEnergyBBDep(const string& line) {
	vector<string> spt;
	splitString(line, " ", &spt);
	this->rotName = spt.at(0);
	for(int i=0;i<200;i++)
	{
		this->energyList[i] = atof(spt.at(i+1).c_str());
	}
}

RotamerEnergyBBDep::~RotamerEnergyBBDep() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPdesignseq */
