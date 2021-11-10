/*
 * ResName.h
 *
 *  Created on: 2017��10��17��
 *      Author: notxp
 */

#ifndef DESIGNSEQ_RESNAME_H_
#define DESIGNSEQ_RESNAME_H_

#include <string>
#include <vector>
#include <map>

namespace NSPdesignseq {

using namespace std;
class ResName {
private:
	string aaSeq;
	vector<string> triNameList;
	map<string,int> triToIntMap;
	map<char,int> sinToIntMap;

public:
	ResName();
	static const ResName &resname();
	char intToSin(int i) const;
	int sinToInt(char sin) const;
	char triToSin(const string& tri) const;
	string sinToTri(char sin) const;
	int triToInt(const string& tri) const;
	string intToTri(int i) const;
	bool isStandardAminoAcid(const string& tri) const;
	bool isAminoAcid(const string& tri) const;
	int chiNum(const string& tri) const;
	virtual ~ResName();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_RESNAME_H_ */
