/*
 * DesignTemplate.h
 *
 *  Created on: 2017Äê12ÔÂ2ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_DESIGNTEMPLATE_H_
#define DESIGNSEQ_DESIGNTEMPLATE_H_
#include <time.h>
#include "designseq/ProteinRep.h"
#include "designseq/Conformer.h"
#include "designseq/S1EnergyTable.h"
#include "designseq/S2MatrixFinder.h"
#include "designseq/SeqProfile.h"
#include "designseq/Rotamer.h"
#include "designseq/RotamerLib.h"
#include "designseq/EnergyArray.h"
#include "designseq/EnergyMatrix.h"
#include "designseq/DesignParameters.h"
#include "designseq/StructureInfo.h"
#include "designseq/AtomicEnergyCalcular.h"
#include "designseq/AtomLib.h"
#include "designseq/AADesignTemplate.h"
#include "backbone/backbonesite.h"

namespace NSPdesignseq {

using namespace NSPproteinrep;

class EnergyCalculatorTemplate{
public:
	vector<EnergyArray*> eaList;
	vector<EnergyMatrix*> emList;
	vector<BackboneSitesPair*> bsPairs;
	vector<vector<int>*> involvedPairs;

};

class DesignTemplate: public EnergyCalculatorTemplate {

public:
//	ProteinChain* pc;
	int resNum;
//	vector<Residue*> resList;

	map<string, int> chainIDResIDToIndex;

	vector<BackBoneSite*> bsList;
//	vector<BackboneSitesPair*> bsPairs;

//	vector<vector<int>*> involvedPairs;
	vector<vector<int>*> neighborResidues;

	vector<AAProbabilityArray*> paList;
	vector<AAScoreMatrix*> smList;

	RotamerLib* rotLibA;
	RotamerLib* rotLibB;

	AtomLib atLib;

	vector<RotamerGroup*> rotGroups;

//	vector<EnergyArray*> eaList;
//	vector<EnergyMatrix*> emList;

	AtomicEnergyCalcular aec;
	ResName rn;
	DesignParameters* dp;

	S1EnergyTable* s1ET;



	DesignTemplate(vector<BackBoneSite*>& bsList, DesignParameters* dp, S1EnergyTable& s1Etable, S2MatrixFinder& s2Etable);
	DesignTemplate(vector<BackBoneSite*>& bsList, string& resFile, DesignParameters* dp, S1EnergyTable& s1Etable, S2MatrixFinder& s2Etable);



	void checkResInvolvedMap(){
		for(int i=0;i<this->involvedPairs.size();i++){
			int involvedNum = this->involvedPairs[i]->size();
			printf("pos: %3d pairNum: %3d\n",i,involvedNum);
			for(int j=0;j<involvedNum;j++){
				int pairID = this->involvedPairs[i]->at(j);
				BackboneSitesPair* pair = this->bsPairs.at(pairID);
				int posA = pair->siteA->resseq;
				int posB = pair->siteB->resseq;
				printf("pair: %3d %3d\n",posA, posB);
			}

		}
	}

	void loadSingleResidueEnergy();

	void loadPairwiseEnergy();
	void updateResPairs();
	void loadS1S2(S1EnergyTable& s1Etable, S2MatrixFinder& s2Etable);
	void checkS1(){
		for(int i=0;i<resNum;i++){
			int choiceNum = this->rotGroups.at(i)->rotNum;
			int x = this->eaList.at(i)->getChoiceNum();
			cout << "I: " << i << " " << choiceNum << " " << x << endl;
		}
	}
	void getSeqProfile(S1EnergyTable& s1Etable, S2MatrixFinder& s2Etable, SeqProfile* result);


	bool contact(int posA, int posB, float cutoff);
	void addProfileRotamerGroups(SeqProfile* prof, float cutoff);
	void addRotamerGroupsFromResFile(string& resFile, SeqProfile* prof, float cutoff);
	float clashScore(vector<Atom*>& listA, vector<Atom*>& listB, string typeB);
	void deleteSc_BB_ClashedRotamers(float cutoff);
	void printRotChoice();
	virtual ~DesignTemplate();
};

float pairEnergy(Conformer* confA, Conformer* confB,int seqSep, float wt, AtomicEnergyCalcular* aec);

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_DESIGNTEMPLATE_H_ */
