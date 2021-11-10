/*
 * AASequence.h
 *
 *  Created on: 2017Äê12ÔÂ15ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_AASEQUENCE_H_
#define DESIGNSEQ_AASEQUENCE_H_
#include <string>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include "designseq/ResName.h"
namespace NSPdesignseq {
using namespace std;

class AASequence {
private:
	int len;
	int* seqInt;
public:
	AASequence(string seq);
	AASequence(int* seqInts, int len);
	AASequence(int len);
	void copyChoice(AASequence* clone);
	void getMutSeq(int pos, int newChoice, AASequence* mutSeq);
	void applyMutation(int pos, int newChoice);
	int getChoice(int pos);
	int getLength();
	string toString();
	virtual ~AASequence();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_AASEQUENCE_H_ */
