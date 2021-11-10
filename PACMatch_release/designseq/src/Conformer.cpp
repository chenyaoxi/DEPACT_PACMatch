/*
 * Conformer.cpp
 *
 *  Created on: 2018Äê1ÔÂ10ÈÕ
 *      Author: notxp
 */

#include "designseq/Conformer.h"

namespace NSPdesignseq {

Conformer::~Conformer(){
	XYZ* t;
	for(int i=0;i<this->bbCoordList.size();i++){
		t = this->bbCoordList.at(i);
		delete t;
	}
	for(int i=0;i<this->scCoordList.size();i++){
		t = this->scCoordList.at(i);
		delete t;
	}

}


ConformerGroup::~ConformerGroup(){
	for(int i=0;i<this->confList.size();i++){
		delete confList.at(i);
	}
}
} /* namespace NSPdesignseq */
