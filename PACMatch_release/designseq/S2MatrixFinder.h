/*
 * S2MatrixFinder.h
 *
 *  Created on: 2017Äê10ÔÂ26ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_S2MATRIXFINDER_H_
#define DESIGNSEQ_S2MATRIXFINDER_H_

#include "designseq/ProteinRep.h"
#include "designseq/AAScoreMatrix.h"
#include <map>
#include <vector>
#include <string>

namespace NSPdesignseq {
using namespace std;

class S2MatrixFinder {
private:
	map<string,vector<ResPairOrientation>> repRPs;
	map<string,vector<SaiPair>> repSAIs;

	void loadRepPoints(string key, string fileName);
	void loadAllRepPoints(string& path);
	void loadSaiPoints(string key, string fileName);
	void loadAllSaiPoints(string& path);
	string findMatrixFile(BackboneSitesPair& bsPair);
	int findSaiIndex(BackboneSitesPair& bsPair);
	int findRPOIndex(BackboneSitesPair& bsPair);


public:
	S2MatrixFinder();

	void getSM(BackboneSitesPair* rp, AAScoreMatrix* outputSM);
	virtual ~S2MatrixFinder();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_S2MATRIXFINDER_H_ */
