/*
 * AADesignTemplate.h
 *
 *  Created on: 2017Äê12ÔÂ15ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_AADESIGNTEMPLATE_H_
#define DESIGNSEQ_AADESIGNTEMPLATE_H_

#include "designseq/ProteinRep.h"
#include "designseq/S1EnergyTable.h"
#include "designseq/S2MatrixFinder.h"
#include "designseq/Rotamer.h"
#include "designseq/RotamerLib.h"
#include "designseq/EnergyArray.h"
#include "designseq/EnergyMatrix.h"
#include "designseq/DesignParameters.h"
#include "designseq/StructureInfo.h"
#include "designseq/AtomicEnergyCalcular.h"
#include "designseq/AASequence.h"
#include "designseq/SeqProfile.h"
#include "backbone/backbonesite.h"

namespace NSPdesignseq {
using namespace NSPproteinrep;

class AADesignTemplate {
private:
	vector<BackBoneSite*> bsList;
	vector<BackboneSitesPair*> bsPairs;
	vector<AAProbabilityArray*> paList;
	vector<AAScoreMatrix*> smList;
	vector<vector<int>*> involvedPairs;

public:
	AADesignTemplate(vector<BackBoneSite*>& bsList, S1EnergyTable* s1Etable, S2MatrixFinder* s2Etable);

	float energy(AASequence* seq);
	float mutEnergy(AASequence* oldSeq, int pos, int newChoice);
	bool accept(float deltaE, float T);
	string mcMinimize();
	void fixTemperatureDesign(float T, int seqNumber, vector<string>& results);

	void getProfile(SeqProfile* prof);
	virtual ~AADesignTemplate();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_AADESIGNTEMPLATE_H_ */
