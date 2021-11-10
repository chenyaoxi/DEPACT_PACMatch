/*
 * RotamerEnergyBBDep.h
 *
 *  Created on: 2017Äê10ÔÂ23ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_ROTAMERENERGYBBDEP_H_
#define DESIGNSEQ_ROTAMERENERGYBBDEP_H_

#include <string>
#include <stdlib.h>
#include <iostream>
#include "designseq/StringTool.h"

namespace NSPdesignseq {

using namespace std;

class RotamerEnergyBBDep {
public:
	float energyList[200];
	string rotName;

	RotamerEnergyBBDep(const string& line);
	virtual ~RotamerEnergyBBDep();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_ROTAMERENERGYBBDEP_H_ */
