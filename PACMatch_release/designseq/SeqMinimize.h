/*
 * SeqMinimize.h
 *
 *  Created on: 2017Äê12ÔÂ26ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_SEQMINIMIZE_H_
#define DESIGNSEQ_SEQMINIMIZE_H_

#include <string>
#include <vector>
#include "stdio.h"
#include <iostream>
#include "time.h"
#include "designseq/RotamerGroup.h"
#include "designseq/DesignTemplate.h"

namespace NSPdesignseq {

using namespace std;

class RotSequence {
private:
	int seqLength;
	int* seqChoices;
	vector<RotamerGroup*> rotGroups;

public:
	RotSequence(int len);
	RotSequence(vector<RotamerGroup*>& groups);
	void copyValue(RotSequence* from);
	void setRandomChoice();
	void applyMutation(int pos, int choice);

	string toAASequence();
	int getChoice(int pos) { return this->seqChoices[pos];}
	Rotamer* getRotamer(int pos, int choice) {return this->rotGroups[pos]->rotList.at(choice);}
	Rotamer* getRotamer(int pos) {
		int choice = this->seqChoices[pos];
		return this->rotGroups[pos]->rotList.at(choice);
	}
	int getLength() {return this->seqLength;}
	int choiceNum(int pos) {return this->rotGroups[pos]->rotNum;}

	RotamerGroup* getRotGroup(int pos) {return this->rotGroups[pos];}

	double mutEnergy(int pos, int mutChoice, EnergyCalculatorTemplate* ec);
	double resInvolvedEnergy(int pos, int choice, EnergyCalculatorTemplate* ec);
	void acceptMut(int pos, int mutChoice);
	double totEnergy(EnergyCalculatorTemplate* ec);
	void printDetailEnergy(EnergyCalculatorTemplate* ec);

	void checkChoice(){
		for(int i=0;i<seqLength;i++){
			RotamerGroup* group = rotGroups.at(i);
			int rotNum = group->rotNum;
			int choice = this->seqChoices[i];
			cout << "site: " << i << " " << rotNum << " " << choice << endl;
		}
	}
	virtual ~RotSequence();
};



class DesignMC {
private:
	DesignTemplate* dt;
	double T0;
	double T1;
	double annealFactor;
	int step;

public:
	DesignMC(DesignTemplate* dt){
		this->dt = dt;
		this->T0 = 200.0;
		this->T1 = 0.01;
		this->annealFactor = 0.95;
		this->step = 50000;
		srand((unsigned)time(NULL));
	}

	bool accept(double mutEnergy, double T);

	void mcRun(RotSequence* result);

	void printPDB(RotSequence* unit, string outputFile);

	virtual ~DesignMC();
};



} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_SEQMINIMIZE_H_ */
