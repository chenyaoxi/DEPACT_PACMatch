/*
 * S1EnergyTable.h
 *
 *  Created on: 2017Äê10ÔÂ26ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_S1ENERGYTABLE_H_
#define DESIGNSEQ_S1ENERGYTABLE_H_

#include "designseq/ProteinRep.h"
#include "backbone/backbonesite.h"
#include "dataio/datapaths.h"
#include "designseq/StringTool.h"
#include "designseq/AAProbabilityArray.h"

namespace NSPdesignseq {

class S1EnergyTable {

private:
	float s1ETable[3][18][200][20];
	PhipsiLib ppLib;

	int ssToInt(char c){
		if(c == 'H' || c == 'G' || c == 'I') return 0;
		else if(c == 'E') return 1;
		else return 2;
	}



public:
	S1EnergyTable();
	int saiToInt(float sai){
		if(sai < 0.15) return 0;
		else if(sai < 0.25) return 1;
		else if(sai < 0.28) return 2;
		else return (int)(sai*20-2.5);
	}
	void getS1(NSPproteinrep::BackBoneSite& bbSite, AAProbabilityArray* pa);
	virtual ~S1EnergyTable();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_S1ENERGYTABLE_H_ */
