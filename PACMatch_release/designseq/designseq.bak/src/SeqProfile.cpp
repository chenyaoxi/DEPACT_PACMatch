/*
 * SeqProfile.cpp
 *
 *  Created on: 2017Äê12ÔÂ16ÈÕ
 *      Author: notxp
 */

#include "designseq/SeqProfile.h"

namespace NSPdesignseq {

SeqProfile::SeqProfile(vector<string>& seqs) {
	initProfile(seqs);
}

void SeqProfile::initProfile(vector<string>& seqs){

	if(seqs.size() == 0){
		cout << "no sequences in profile" << endl;
		exit(1);
	}

	this->seqLen = seqs.at(0).length();
	for(string s : seqs){
		this->seqs.push_back(s);
	}

	int counts[seqLen][20];
	for(int i=0;i<seqLen;i++){
		for(int j=0;j<20;j++)
			counts[i][j] = 0;
	}

	ResName rn;
	for(int i=0;i<seqs.size();i++){
		string s = seqs.at(i);
		//cout << "i: " << i << " " << s << endl;
		for(int j=0;j<seqLen;j++){
			char c = s.at(j);
			int aa = rn.sinToInt(c);
			counts[j][aa] ++;
		}
	}

	for(int i=0;i<seqLen;i++){
		AAProbabilityArray pa(counts[i]);
		this->paList.push_back(pa);
	}

}

float SeqProfile::getP(int pos, int aa){
	return this->paList.at(pos).getProbability(aa);
}

float SeqProfile::getRelP(int pos, int aa){
	float bg[20] = 	{0.0839, 0.0126, 0.0593, 0.0681, 0.0409,
            0.0735, 0.0235, 0.0576, 0.0567, 0.0940,
            0.0167, 0.0427, 0.0463, 0.0375, 0.0517,
            0.0588, 0.0546, 0.0713, 0.0146, 0.0358};
	float p = this->paList.at(pos).getProbability(aa);
	return p/bg[aa];
}

float SeqProfile::getProfScore(int pos, int aa){
	float bg[20] = 	{0.0839, 0.0126, 0.0593, 0.0681, 0.0409,
	            0.0735, 0.0235, 0.0576, 0.0567, 0.0940,
	            0.0167, 0.0427, 0.0463, 0.0375, 0.0517,
	            0.0588, 0.0546, 0.0713, 0.0146, 0.0358};
	float p = this->paList.at(pos).getProbability(aa);
	float score = 5.0*log(p/bg[aa]);
	return score;
}

void SeqProfile::printProfile(){

	this->paList.at(0).printArray();

	cout << "    ";
	ResName rn;
	for(int i=0;i<20;i++){
		cout << "  " << rn.intToTri(i) << " ";
	}
	cout << endl;
	for(int i=0;i<seqLen;i++){
		printf("%-3d ",i);
		for(int j=0;j<20;j++){
			printf(" %5.3f",getP(i,j));
		}
		cout << endl;
	}
}

SeqProfile::~SeqProfile() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPdesignseq */
