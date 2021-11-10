/*
 * EnergyArray.cpp
 *
 *  Created on: 2017Äê10ÔÂ17ÈÕ
 *      Author: notxp
 */

#include "designseq/EnergyArray.h"

namespace NSPdesignseq {

EnergyArray::EnergyArray(int choiceNum) {
	this->ea = new float[choiceNum];
	this->choiceNum = choiceNum;
}

int EnergyArray::getChoiceNum() {
	return this->choiceNum;
}

void EnergyArray::setEnergy(int choice, float e) {
	if(choice < 0 || choice>=choiceNum)
	{
		cerr << "choiceNum: " << choiceNum << " invalid choice: " << choice << endl;
		exit(0);
	}
	this->ea[choice] = e;
}

float EnergyArray::getEnergy(int choice){
	if(choice < 0 || choice>=choiceNum)
	{
		cerr << "choiceNum: " << choiceNum << " invalid choice: " << choice << endl;
		exit(0);
	}
	return this->ea[choice];
}

EnergyArray::~EnergyArray() {
	delete [] ea;
}

} /* namespace NSPdesignseq */
