/*
 * EnergyMatrix.cpp
 *
 *  Created on: 2017Äê10ÔÂ17ÈÕ
 *      Author: notxp
 */

#include "designseq/EnergyMatrix.h"

namespace NSPdesignseq {

EnergyMatrix::EnergyMatrix(int choiceANum, int choiceBNum, int posA, int posB) {
	this->choiceANum = choiceANum;
	this->choiceBNum = choiceBNum;
	this->posA = posA;
	this->posB = posB;
	this->em = new float[choiceANum*choiceBNum];
	int n=choiceANum*choiceBNum;
	for(int i=0;i<n;i++)
	{
		this->em[i] = 0;
	}
}

EnergyMatrix::~EnergyMatrix() {
	// TODO Auto-generated destructor stub
	delete [] em;
}

} /* namespace NSPdesignseq */
