/*
 * RotamerGroup.cpp
 *
 *  Created on: 2017Äê10ÔÂ23ÈÕ
 *      Author: notxp
 */

#include "designseq/RotamerGroup.h"

namespace NSPdesignseq {

RotamerGroup::RotamerGroup() {
	// TODO Auto-generated constructor stub
	for(int i=0;i<20;i++){
		aaRots.push_back(new vector<Rotamer*>());
	}
	this->rotNum = 0;
}


void RotamerGroup::addRotamer(Rotamer* rot)
{
	ResName rn;
	int aaType = rn.triToInt(rot->triName);
	this->rotList.push_back(rot);
	this->aaRots.at(aaType)->push_back(rot);
	this->rotNum ++;
}

void RotamerGroup::deleteRotamer(int i){
	Rotamer* rot = rotList.at(i);
	string rotName = rot->rotName;
	ResName rn;
	int aaType = rn.triToInt(rot->triName);
	vector<Rotamer*>* aaRotList = aaRots.at(aaType);
	for(int k=0;k<aaRotList->size();k++){
		Rotamer* r = aaRotList->at(k);
		if(r->rotName == rot->rotName){
			aaRotList->erase(aaRotList->begin()+k);
			break;
		}
	}
	this->rotList.erase(rotList.begin()+i);
	this->rotNum --;
}

void RotamerGroup::clear(){
	this->rotList.clear();
	this->rotNum = 0;
}

RotamerGroup::~RotamerGroup(){
	for(int i=0;i<20;i++){
		delete aaRots.at(i);
	}
}


} /* namespace NSPdesignseq */
