/*
 * S2MatrixFinder.cpp
 *
 *  Created on: 2017Äê10ÔÂ26ÈÕ
 *      Author: notxp
 */

#include "designseq/S2MatrixFinder.h"

namespace NSPdesignseq {

S2MatrixFinder::S2MatrixFinder() {
	string repPath = NSPdataio::datapath() + "S2EnergyTable/repPoints/";

	loadAllRepPoints(repPath);

	loadAllSaiPoints(repPath);

}

void S2MatrixFinder::loadRepPoints(string key, string fileName){

     ifstream file;
     file.open(fileName.c_str(), ios::in);
     if(! file.is_open()){
    	 cout << "fail to open file " << fileName << endl;
    	 exit(1);
     }
     string s;
     string substr;
     vector<ResPairOrientation> rpList;
     while(getline(file,s)){
    	 if(s.length() < 100) continue;
    	 substr = s.substr(8,119);
    	 ResPairOrientation rpo(substr);
    	 rpList.push_back(rpo);
     }
     repRPs[key] = rpList;
}


void S2MatrixFinder::loadAllRepPoints(string& path){
	loadRepPoints("HH1", path + "HH1-100");
	loadRepPoints("HE1", path + "HE1-100");
	loadRepPoints("HC1", path + "HC1-200");
	loadRepPoints("EH1", path + "EH1-100");
	loadRepPoints("EE1", path + "EE1-100");
	loadRepPoints("EC1", path + "EC1-200");
	loadRepPoints("CH1", path + "CH1-200");
	loadRepPoints("CE1", path + "CE1-200");
	loadRepPoints("CC1", path + "CC1-200");

	loadRepPoints("HH2", path + "HH2-100");
	loadRepPoints("HE2", path + "HE2-100");
	loadRepPoints("HC2", path + "HC2-300");
	loadRepPoints("EH2", path + "EH2-100");
	loadRepPoints("EE2", path + "EE2-100");
	loadRepPoints("EC2", path + "EC2-400");
	loadRepPoints("CH2", path + "CH2-300");
	loadRepPoints("CE2", path + "CE2-400");
	loadRepPoints("CC2", path + "CC2-400");

	loadRepPoints("HH3", path + "HH3-100");
	loadRepPoints("HE3", path + "HE3-100");
	loadRepPoints("HC3", path + "HC3-400");
	loadRepPoints("EH3", path + "EH3-100");
	loadRepPoints("EE3", path + "EE3-200");
	loadRepPoints("EC3", path + "EC3-600");
	loadRepPoints("CH3", path + "CH3-400");
	loadRepPoints("CE3", path + "CE3-600");
	loadRepPoints("CC3", path + "CC3-600");

	loadRepPoints("HH4", path + "HH4-100");
	loadRepPoints("HE4", path + "HE4-200");
	loadRepPoints("HC4", path + "HC4-800");
	loadRepPoints("EH4", path + "EH4-200");
	loadRepPoints("EE4", path + "EE4-200");
	loadRepPoints("EC4", path + "EC4-800");
	loadRepPoints("CH4", path + "CH4-800");
	loadRepPoints("CE4", path + "CE4-800");
	loadRepPoints("CC4", path + "CC4-1000");

	loadRepPoints("HH5", path + "HH5-2000");
	loadRepPoints("HE5", path + "HE5-2000");
	loadRepPoints("HC5", path + "HC5-3000");
	loadRepPoints("EH5", path + "EH5-2000");
	loadRepPoints("EE5", path + "EE5-2000");
	loadRepPoints("EC5", path + "EC5-3000");
	loadRepPoints("CH5", path + "CH5-3000");
	loadRepPoints("CE5", path + "CE5-3000");
	loadRepPoints("CC5", path + "CC5-3000");
}

void S2MatrixFinder::loadSaiPoints(string key, string fileName){
	vector<SaiPair> saiList;
    ifstream file;
    file.open(fileName.c_str(), ios::in);
    if(! file.is_open()){
   	 cout << "fail to open file " << fileName << endl;
   	 exit(1);
    }
    string s;


    vector<string> spt;
    while(getline(file,s)){
   	 if(s.length() < 10) continue;
   	 splitString(s, " ", &spt);
   	 double x = atof(spt[0].c_str());
   	 double y = atof(spt[1].c_str());
   	 SaiPair sp(x,y);
   	 saiList.push_back(sp);
    }
    this->repSAIs[key] = saiList;
}

void S2MatrixFinder::loadAllSaiPoints(string& path){



	string hec = "HEC";
	char tmp[100];
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			for(int k=1;k<5;k++)
			{
				sprintf(tmp, "%c%c%d", hec.at(i), hec.at(j), k);
				string key = string(tmp);
				sprintf(tmp, "%spair%d-sai",path.c_str(),k);
				string fileName = string(tmp);
				loadSaiPoints(key, fileName);
			}


			sprintf(tmp, "%c%c%d", hec.at(i), hec.at(j), 5);
			string key = string(tmp);
			string fileName = path + key + "-sai";
			loadSaiPoints(key, fileName);
		}
	}
}

string S2MatrixFinder::findMatrixFile(BackboneSitesPair& bsPair){
	char ssA = bsPair.siteA->sscode;
	char ssB = bsPair.siteB->sscode;
	int saiIndex = findSaiIndex(bsPair) + 1;
	char tmp[100];
	sprintf(tmp, "S2-%c%c%d-%d-1500", ssA, ssB, bsPair.seqSep, saiIndex);
	string key = string(tmp);

	return NSPdataio::datapath() + "S2EnergyTable/strict/" + key;
}

int S2MatrixFinder::findSaiIndex(BackboneSitesPair& bsPair){
//	cout << "find sai index: " << endl;
	char ssA = bsPair.siteA->sscode;
	char ssB = bsPair.siteB->sscode;
	int seqSep = bsPair.seqSep;
	char tmp[20];
	sprintf(tmp, "%c%c%d", ssA, ssB, seqSep);
	string key = string(tmp);
	double saiA = bsPair.siteA->data_[3];
	double saiB = bsPair.siteB->data_[3];
	SaiPair sp(saiA, saiB);

//	cout << "sai key: " << key << endl;


	vector<SaiPair>& list = this->repSAIs[key];
//	cout << "list size: " << list.size() << endl;

	double minDist = 10;
	int minIndex = -1;
	for(int i=0;i<12;i++){
		double d = sp.distanceSquare(&(list.at(i)));
		if(d < minDist){
			minDist = d;
			minIndex = i;
		}
	}
	return minIndex;


}

int S2MatrixFinder::findRPOIndex(BackboneSitesPair& bsPair){
	char ssA = bsPair.siteA->sscode;
	char ssB = bsPair.siteB->sscode;
	int seqSep = bsPair.seqSep;
	char tmp[20];
	sprintf(tmp, "%c%c%d", ssA, ssB, seqSep);
	string key = string(tmp);

	ResPairOrientation rpo = bsPair.ori;


	vector<ResPairOrientation>& list = this->repRPs[key];





	double minDist = 10000;
	int minIndex = -1;
	for(int i=0;i<list.size();i++){
		double d = rpo.rmsd(list.at(i));
		if(d < minDist){
			minDist = d;
			minIndex = i;
		}
	}
	return minIndex;
}

void S2MatrixFinder::getSM(BackboneSitesPair* rp, AAScoreMatrix* outputSM){

//	cout << "begin" << endl;
	string fileName = findMatrixFile(*rp);
//	cout << "fileName: " << fileName << endl;

	int rpoIndex = findRPOIndex(*rp);
	if(rpoIndex < 0){
		cerr << "can't find rp " << endl;
		exit(1);
	}

//	cout << "rpoIndex: " << rpoIndex << endl;


	int start = rpoIndex*22+2;
//	cout << "start : " << start << endl;
	ifstream file;
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()){
	       cout << "fail to open file " << fileName << endl;
	       exit(1);
	}
	string s;

	for(int i=0;i<start;i++)
		getline(file, s);
	vector<string> spt;
	float matrix[20][20];
	for(int i=0;i<20;i++){
		getline(file,s);
		splitString(s," ",&spt);
		for(int k=0;k<20;k++){
			matrix[i][k] = atof(spt[k].c_str());
		}
	}
	file.close();
	outputSM->initValue(matrix);

}

S2MatrixFinder::~S2MatrixFinder() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPdesignseq */
