/*
 * AASequence.cpp
 *
 *  Created on: 2017��12��15��
 *      Author: notxp
 */

#include "designseq/AASequence.h"

namespace NSPdesignseq {
AASequence::AASequence(string seq){
	ResName rn;
	this->len = seq.length();
	this->seqInt = new int[len];
	for(int i=0;i<this->len;i++){
		seqInt[i] = rn.sinToInt(seq.at(i));
	}
}

AASequence::AASequence(int seqInts[], int len){
	this->len = len;
	this->seqInt = new int[len];
	for(int i=0;i<len;i++){
		this->seqInt[i] = seqInts[i];
	}
}

AASequence::AASequence(int len) {
	this->len = len;
	this->seqInt = new int[len];
	srand((unsigned)time(NULL));
	for(int i=0;i<len;i++){
		this->seqInt[i] = rand()%20;
	}
}

void AASequence::copyChoice(AASequence* clone){
	if(this->len != clone->len){
		cout << "sequence length not equal"	<< endl;
		exit(1);
	}
	for(int i=0;i<this->len;i++){
		clone->seqInt[i] = this->seqInt[i];
	}
}

void AASequence::getMutSeq(int pos, int newChoice, AASequence* mutSeq){
	if(this->len != mutSeq->len){
		cout << "sequence length not equal"	<< endl;
		exit(1);
	}
	for(int i=0;i<this->len;i++){
		mutSeq->seqInt[i] = this->seqInt[i];
	}
	mutSeq->seqInt[pos] = newChoice;
}

void AASequence::applyMutation(int pos, int newChoice){
	this->seqInt[pos] = newChoice;
}

int AASequence::getChoice(int pos){
	return this->seqInt[pos];
}

int AASequence::getLength(){
	return this->len;
}

string AASequence::toString(){
	char s[this->len+1];
	ResName rn;
	for(int i=0;i<this->len;i++){
		s[i] = rn.intToSin(this->seqInt[i]);
	}
	s[this->len]='\0';
	string seq = s;
	return seq;
}
AASequence::~AASequence() {
	// TODO Auto-generated destructor stub
	delete []seqInt;
}

} /* namespace NSPdesignseq */
