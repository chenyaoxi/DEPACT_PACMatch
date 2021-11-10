/*
 * EnergyMatrix.h
 *
 *  Created on: 2017Äê10ÔÂ17ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_ENERGYMATRIX_H_
#define DESIGNSEQ_ENERGYMATRIX_H_

#include <iostream>
#include "stdlib.h"

namespace NSPdesignseq {

using namespace std;

class EnergyMatrix {
public:
	int choiceANum;
	int choiceBNum;
	int posA, posB;
	float* em;
	EnergyMatrix(int choiceANum, int choiceBNum, int posA, int posB);
	void setEnergy(int choiceA, int choiceB, float e){
		if(choiceA <0 || choiceA >= this->choiceANum)
		{
			cerr << "choiceA out of boundary: " << choiceA << endl;
			exit(1);
		}

		if(choiceB <0 || choiceB >= this->choiceBNum)
		{
			cerr << "choiceB out of boundary: " << choiceB << endl;
			exit(1);
		}
		this->em[choiceA*choiceBNum + choiceB] = e;
	};

	float getEnergy(int choiceA, int choiceB)
	{
		if(choiceA <0 || choiceA >= this->choiceANum)
		{
			cerr << "POSA: " << posA << " posB " << posB << endl;
			cerr << "EM GET_ENERGY: choiceA out of boundary: " << choiceA << " pos: " << posA <<  " totChoiceNum: " << this->choiceANum <<  endl;
			exit(1);
		}

		if(choiceB <0 || choiceB >= this->choiceBNum)
		{
			cerr << "EM GET_ENERGY: choiceB out of boundary: " << choiceB << " choiceNum: " << choiceBNum << endl;
			exit(1);
		}

		int n=choiceANum*choiceBNum;
		int m = choiceA*choiceBNum + choiceB;
		if(m >= n || m < 0) {
			cerr << "choiceA: " << choiceA << " choiceB: " << choiceB << " n: " << n << " m: " << m << endl;
			abort();
		}
		return this->em[choiceA*choiceBNum + choiceB];
	}
	virtual ~EnergyMatrix();
	};



}
		 /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_ENERGYMATRIX_H_ */
