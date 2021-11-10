/*
 * ResName.cpp
 *
 *  Created on: 2017��10��17��
 *      Author: notxp
 */

#include "designseq/ResName.h"

namespace NSPdesignseq {
const ResName &ResName::resname(){
	static ResName stdname;
	return stdname;
}
ResName::ResName() {
	// TODO Auto-generated constructor stub
	aaSeq = "ACDEFGHIKLMNPQRSTVWYBJX";

	this->triNameList.reserve(23);
	this->triNameList.push_back("ALA");
	this->triNameList.push_back("CYS");
	this->triNameList.push_back("ASP");
	this->triNameList.push_back("GLU");
	this->triNameList.push_back("PHE");
	this->triNameList.push_back("GLY");
	this->triNameList.push_back("HIS");
	this->triNameList.push_back("ILE");
	this->triNameList.push_back("LYS");
	this->triNameList.push_back("LEU");
	this->triNameList.push_back("MET");
	this->triNameList.push_back("ASN");
	this->triNameList.push_back("PRO");
	this->triNameList.push_back("GLN");
	this->triNameList.push_back("ARG");
	this->triNameList.push_back("SER");
	this->triNameList.push_back("THR");
	this->triNameList.push_back("VAL");
	this->triNameList.push_back("TRP");
	this->triNameList.push_back("TYR");
	this->triNameList.push_back("MSE");
	this->triNameList.push_back("PSD");
	this->triNameList.push_back("UNK");

	this->triToIntMap["ALA"] = 0;
	this->triToIntMap["CYS"] = 1;
	this->triToIntMap["ASP"] = 2;
	this->triToIntMap["GLU"] = 3;
	this->triToIntMap["PHE"] = 4;
	this->triToIntMap["GLY"] = 5;
	this->triToIntMap["HIS"] = 6;
	this->triToIntMap["ILE"] = 7;
	this->triToIntMap["LYS"] = 8;
	this->triToIntMap["LEU"] = 9;
	this->triToIntMap["MET"] = 10;
	this->triToIntMap["ASN"] = 11;
	this->triToIntMap["PRO"] = 12;
	this->triToIntMap["GLN"] = 13;
	this->triToIntMap["ARG"] = 14;
	this->triToIntMap["SER"] = 15;
	this->triToIntMap["THR"] = 16;
	this->triToIntMap["VAL"] = 17;
	this->triToIntMap["TRP"] = 18;
	this->triToIntMap["TYR"] = 19;
	this->triToIntMap["MSE"] = 20;
	this->triToIntMap["PSD"] = 21;
	this->triToIntMap["UNK"] = 22;

	this->sinToIntMap['A'] = 0;
	this->sinToIntMap['C'] = 1;
	this->sinToIntMap['D'] = 2;
	this->sinToIntMap['E'] = 3;
	this->sinToIntMap['F'] = 4;
	this->sinToIntMap['G'] = 5;
	this->sinToIntMap['H'] = 6;
	this->sinToIntMap['I'] = 7;
	this->sinToIntMap['K'] = 8;
	this->sinToIntMap['L'] = 9;
	this->sinToIntMap['M'] = 10;
	this->sinToIntMap['N'] = 11;
	this->sinToIntMap['P'] = 12;
	this->sinToIntMap['Q'] = 13;
	this->sinToIntMap['R'] = 14;
	this->sinToIntMap['S'] = 15;
	this->sinToIntMap['T'] = 16;
	this->sinToIntMap['V'] = 17;
	this->sinToIntMap['W'] = 18;
	this->sinToIntMap['Y'] = 19;
	this->sinToIntMap['B'] = 20;
	this->sinToIntMap['J'] = 21;
	this->sinToIntMap['X'] = 22;
}

char ResName::intToSin(int i) const{
	if(i>=0 && i<=22)
			return aaSeq.at(i);
		return 'X';
}

int ResName::sinToInt(char sin) const{
	map<char,int>::const_iterator it = sinToIntMap.find(sin);
	if(it != sinToIntMap.end())
		return it->second;
	else
		return 22;
}

int ResName::triToInt(const string& tri) const{
	map<string,int>::const_iterator it = triToIntMap.find(tri);
	if(it != triToIntMap.end())
		return it->second;

	if(tri == "HID" || tri == "HIE" || tri == "HIP")
		return 6;
	return 22;
}

string ResName::intToTri(int i) const{
	if(i>=0 && i<=22)
		return triNameList.at(i);
	else
		return "UNK";
}

char ResName::triToSin(const string& tri) const{
	return intToSin(triToInt(tri));
}

string ResName::sinToTri(char sin) const{
	return intToTri(sinToInt(sin));
}

bool ResName::isStandardAminoAcid(const string& tri) const {
	return triToInt(tri) < 20;
}

bool ResName::isAminoAcid(const string& tri) const{
    int aaInt = triToInt(tri);
    if(aaInt < 22) return true;
    else if(tri == "GLX" || tri == "LEF" || tri == "PCA" || tri == "OCS") return true;
    else return false;

}

int ResName::chiNum(const string& tri) const
{
	int chiNum[] = {0, 1, 2, 3, 2, 0, 2, 2, 4, 2, 3, 2, 1, 3, 4, 1, 1, 1, 2, 2};
	int type = triToInt(tri);
	if(type > 19)
		return 0;
	else
		return chiNum[type];

}


ResName::~ResName() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPdesignseq */
