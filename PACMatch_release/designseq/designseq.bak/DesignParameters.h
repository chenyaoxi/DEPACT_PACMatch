/*
 * DesignParameters.h
 *
 *  Created on: 2017Äê12ÔÂ2ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_DESIGNPARAMETERS_H_
#define DESIGNSEQ_DESIGNPARAMETERS_H_

#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "designseq/StringTool.h"
#include "dataio/datapaths.h"
#include "designseq/ResName.h"

namespace NSPdesignseq {
using namespace std;

class DesignParameters {
public:

	double wtS1;
	double wtS2Local;
	double wtS2Nonlocal;
	double wtS2EE;

	double wtCorePacking;
	double wtSurfPacking;
	double packSigma;
	double packSAI0;



	double nonHBondPolarWeight;
	double repWeight;

	double ref[20];

	double wtRotEnergy;

	double vdwLamda1;
	double vdwLamda2;

	double wdNP_NP;
	double wdP_NP;
	double wdP_P;

	double wdHB;
	double hbAngleFactor;


	double aroWD;

	DesignParameters();
	DesignParameters(string& paraFile);
	void readRef(string& refFile){
		ifstream f;
		f.open(refFile.c_str(),ios::in);
		string s;
		vector<string> spt;
		ResName rn;
		while(getline(f,s)){
			splitString(s," ",&spt);
			int aa = rn.triToInt(spt[0]);
			float e = atof(spt[1].c_str());
			this->ref[aa] = e;
		}
		f.close();
	}
	double getVdwWeight(double sai){
		return wtSurfPacking + ( wtCorePacking - wtSurfPacking)/(1+exp((sai-packSAI0)/packSigma));
	}
	virtual ~DesignParameters();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_DESIGNPARAMETERS_H_ */
