/*
 * S1EnergyTable.cpp
 *
 *  Created on: 2017Äê10ÔÂ26ÈÕ
 *      Author: notxp
 */

#include "designseq/S1EnergyTable.h"


namespace NSPdesignseq {



S1EnergyTable::S1EnergyTable() {

	string path = NSPdataio::datapath();
	string s1EnergyFile = path+"S1EnergyTable6";
	ifstream file;
	file.open(s1EnergyFile.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open file " << s1EnergyFile << endl;
		exit(1);
	}
	string s;
	vector<string> spt;
	char ssType;
	double sai, phi, psi;
	int ssInt, saiInt, ppInt;
	getline(file,s);
	while(getline(file,s)){
		if(s.length() < 174) continue;
		vector<string> spt;
		splitString(s," ", &spt);
		ssType = s.at(0);
		sai = atof(spt[1].c_str());
		phi = atof(spt[2].c_str());
		psi = atof(spt[3].c_str());
		ssInt = ssToInt(ssType);
		saiInt = saiToInt(sai);
		Phipsi pp(phi, psi);
		ppInt = ppLib.phipsiToIndex(&pp);
		for(int k=0;k<20;k++){
			float p = atof(spt[k+5].c_str());
			if(p == 0) p = 0.00001;
			s1ETable[ssInt][saiInt][ppInt][k] = p;
		}

	}

	// TODO Auto-generated constructor stub
}

void S1EnergyTable::getS1(NSPproteinrep::BackBoneSite& bs, AAProbabilityArray* pa){
	int ss = ssToInt(bs.sscode);
	int sai = saiToInt(bs.data_[3]);
	Phipsi pp(bs.data_[0], bs.data_[1]);
	int ppIndex = ppLib.phipsiToIndex(&pp);
	pa->initProbability(s1ETable[ss][sai][ppIndex]);
}

S1EnergyTable::~S1EnergyTable() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPdesignseq */
