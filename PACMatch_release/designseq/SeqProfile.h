/*
 * SeqProfile.h
 *
 *  Created on: 2017Äê12ÔÂ16ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_SEQPROFILE_H_
#define DESIGNSEQ_SEQPROFILE_H_
#include <string>
#include <vector>
#include "designseq/AAProbabilityArray.h"
#include "designseq/ResName.h"

namespace NSPdesignseq {
using namespace std;

class SeqProfile {
private:
	vector<string> seqs;
	vector<AAProbabilityArray> paList;
	int seqLen;
public:
	SeqProfile(){
		this->seqLen = -1;
	}
	SeqProfile(vector<string>& seqs);
	void initProfile(vector<string>& seqs);

	int getSeqLen(){
		return this->seqLen;
	}
	float getP(int pos, int aa);
	float getRelP(int pos, int aa);
	float getProfScore(int pos, int aa);
	void printProfile();


	virtual ~SeqProfile();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_SEQPROFILE_H_ */
