/*
 * RotamerLib.cpp
 *
 *  Created on: 2017��10��23��
 *      Author: notxp
 */

#include "designseq/RotamerLib.h"
#include <memory>
namespace NSPdesignseq {


RotamerLib & RotamerLib::rotamerlib(const string &rotLibType){
	static std::map<std::string,std::shared_ptr<RotamerLib>> stdlibs;
	   if(stdlibs.find(rotLibType)!= stdlibs.end()) return *(stdlibs.at(rotLibType));
	   else {
		   stdlibs[rotLibType]=std::shared_ptr<RotamerLib>(new RotamerLib(rotLibType));
		   return *(stdlibs.at(rotLibType));
	   }
}
RotamerLib::RotamerLib() {
	// TODO Auto-generated constructor stub

}

RotamerLib::RotamerLib(const string& rotLibType)
{

	if(rotLibType != "A0" && rotLibType != "A1" && rotLibType != "A2" && rotLibType != "A3" && rotLibType != "B0" && rotLibType != "B1" && rotLibType != "B2" && rotLibType != "B3")
	{
		cerr << "invalid rotamer library type: " << rotLibType << endl;
		exit(1);
	}
	ifstream file;
	string rotPath = NSPdataio::datapath()+"rotLib/";
	string rotFileName = rotPath + "rot" + rotLibType + "-p";
	string energyFileName = rotPath + "rot" + rotLibType + "-energy";
	file.open(rotFileName.c_str(), ios::in);
	if(!file.is_open())
	{
		cerr << "fail to open rotamer lib file " << rotFileName << endl;
		exit(1);
	}
	string rotName, resType, atomName, s;
	int atomNum,i, aa, rotIndex=0;
	float rotP, x,y,z;

	for(i=0;i<20;i++)
	{
		RotamerGroup* group = new RotamerGroup();
		aaRotGroups.push_back(group);
	}


	Rotamer* rotGly = new Rotamer("GLY-1", "GLY");
	this->aaRotGroups.at(5)->addRotamer(rotGly);
	this->allRots.addRotamer(rotGly);
	this->rotNameToIndex["GLY-1"] = rotIndex;
	rotIndex ++;


	int rotNum[20] = {0,};
	while(file.good())
	{
		s = rotName;
		file >> rotName >> atomNum >> rotP;
		if(s == rotName)
			break;
		resType = rotName.substr(0,3);
		aa = rn.triToInt(resType);
		rotNum[aa] ++;
		Rotamer* rot = new Rotamer(rotName, resType);

		this->aaRotGroups.at(aa)->addRotamer(rot);
		this->allRots.addRotamer(rot);
		this->rotNameToIndex[rotName] = rotIndex;
		rotIndex ++;

		for(int i=0;i<atomNum;i++)
		{
			file >> atomName >> x >> y >> z;

			rot->addAtom(atomName,NSPgeometry::XYZ(x,y,z));
		}

		rot->updateLawOfRing();
	}

	file.close();


	file.open(energyFileName.c_str(), ios::in);
	if(!file.is_open())
	{
		cerr << "fail to open file: " << energyFileName << endl;
		exit(1);
	}



	while(getline(file,s))
	{
		RotamerEnergyBBDep* rotE = new RotamerEnergyBBDep(s);
		rotEnergyMap[rotE->rotName] = rotE;
	}


	file.close();

}


RotamerGroup* RotamerLib::getAAGroup(string& triName)
{
	int i = rn.triToInt(triName);
	return aaRotGroups.at(i);
}

RotamerLib::~RotamerLib() {
	unsigned int i;
	for(i=0;i<20;i++)
	{
		delete aaRotGroups.at(i);
	}

	vector<Rotamer*>::iterator it;
	string rotName;

	for(it=allRots.rotList.begin();it!=allRots.rotList.end();it++)
	{
		Rotamer* rot = *it;
		rotName = rot->rotName;
		delete rot;
		delete rotEnergyMap[rotName];
	}

	//cout << "rotLib deleted " << endl;

}

float RotamerLib::getRotamerEnergy(string& rotName, int ppType) const
{
	if(rotEnergyMap.find(rotName) == rotEnergyMap.end())
	{
		cerr << "rotamer not defined: " << rotName << endl;
		exit(1);
	}
	if(ppType < 0 || ppType > 199)
	{
		cerr << "invalid phi-psi torsion type: " << ppType << endl;
	}
	const RotamerEnergyBBDep* re = rotEnergyMap.at(rotName);
	return re->energyList[ppType];
}


} /* namespace NSPdesignseq */
