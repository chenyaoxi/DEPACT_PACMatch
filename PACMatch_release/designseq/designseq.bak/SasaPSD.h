/*
 * SasaPSD.h
 *
 *  Created on: 2017Äê10ÔÂ23ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_SASAPSD_H_
#define DESIGNSEQ_SASAPSD_H_
#include <vector>
#include "designseq/ProteinRep.h"

namespace NSPdesignseq {



class SasaPSD {
private:
	XYZ points[120];
	float sasaIndex[121];
	float radii;
	XYZ psdList[3];

public:
	SasaPSD();
	int exposeNum(NSPproteinrep::BackBoneSite* resA, vector<NSPproteinrep::BackBoneSite*>& resList);
	int exposeNum(Residue* resA, vector<Residue*>& resList);
	int exposeNumTest(Residue* resA, vector<Residue*>& resList);


	float pointNumToIndex(int n){
		if(n < 0 || n > 120)
		{
			cerr << "invalid expose point num: " + n << endl;
			exit(1);
		}
		return sasaIndex[n];
	}
	void calSasa(vector<NSPproteinrep::BackBoneSite>&);
	virtual ~SasaPSD();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_SASAPSD_H_ */
