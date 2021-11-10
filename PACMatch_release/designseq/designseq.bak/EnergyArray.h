/*
 * EnergyArray.h
 *
 *  Created on: 2017Äê10ÔÂ17ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_ENERGYARRAY_H_
#define DESIGNSEQ_ENERGYARRAY_H_

#include <iostream>

namespace NSPdesignseq {

using namespace std;

class EnergyArray {
private:
	float* ea;
	int choiceNum;

public:
	EnergyArray(int choiceNum);
	int getChoiceNum();
	void setEnergy(int choice, float e);
	float getEnergy(int choice);
	virtual ~EnergyArray();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_ENERGYARRAY_H_ */
