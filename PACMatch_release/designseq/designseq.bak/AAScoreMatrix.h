/*
 * AAScoreMatrix.h
 *
 *  Created on: 2017Äê10ÔÂ27ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_AASCOREMATRIX_H_
#define DESIGNSEQ_AASCOREMATRIX_H_

#include <string>
#include <iostream>

namespace NSPdesignseq {

using namespace std;
class AAScoreMatrix {
private:
	float sm[20][20];
	float sampleNum;
	string tag;

public:
	AAScoreMatrix(){
		int i,j;
		for(i=0;i<20;i++){
			for(j=0;j<20;j++)
				sm[i][j] = 0;
		}
		this->sampleNum = -1;
		this->tag = "NEW";
	}

	void setValue(int i, int j, float value){
		if(i < 0 || i > 19 || j< 0 || j > 19){
			cout << "invalid sm index" << endl;
			abort();
		}
		this->sm[i][j] = value;
	}

	void setTaG(const string& tag){
		this->tag = tag;
	}

	void setSampleNum(float n){
		this->sampleNum = n;
	}

	void initValue(float sm[20][20]){
		for(int i=0;i<20;i++)
		{
			for(int j=0;j<20;j++)
			{
				this->sm[i][j] = sm[i][j];
			}
		}
	}

	string getTag(){
		return this->tag;
	}

	float getSampleNum(){
		return this->sampleNum;
	}

	float getValue(int i, int j){
		if(i < 0 || i > 19 || j< 0 || j > 19){
			cout << "invalid sm index" << endl;
			abort();
		}
		return this->sm[i][j];
	}

	void multiply(float wt){
		for(int i=0;i<20;i++){
			for(int j=0;j<20;j++){
				this->sm[i][j] = this->sm[i][j]*wt;
			}
		}

	}

	virtual ~AAScoreMatrix();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_AASCOREMATRIX_H_ */
