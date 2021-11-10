/*
 * Conformer.h
 *
 *  Created on: 2018Äê1ÔÂ10ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_CONFORMER_H_
#define DESIGNSEQ_CONFORMER_H_
#include <string>
#include <vector>
#include "designseq/Rotamer.h"
#include "designseq/RotamerGroup.h"
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "designseq/AtomLib.h"
#include "designseq/ProteinRep.h"
#include "backbone/backbonesite.h"

namespace NSPdesignseq {
using namespace std;
using namespace NSPproteinrep;
using namespace NSPgeometry;

class Conformer {
public:
	string triName;
	vector<XYZ*> bbCoordList;
	vector<XYZ*> scCoordList;
	vector<AtomProperty*> bbApList;
	vector<AtomProperty*> scApList;
	bool hasAromaticRing;
	NSPgeometry::XYZ normalVectorOfRing;
	Conformer(){hasAromaticRing = false; triName = "UNK";}
	Conformer(Rotamer* rot, BackBoneSite* bs, AtomLib* atLib){
		this->triName = rot->triName;
		int n = rot->coordList.size();
		LocalFrame cs = getBackboneSiteLocalFrame(*bs);
		string s;
		for(int i=0;i<n;i++){
			s = rot->triName + "-" + rot->atomNameList.at(i);
			scApList.push_back(atLib->getAtomProperty(s));
			XYZ tmp = cs.local2globalcrd(rot->coordList.at(i));
			scCoordList.push_back(new XYZ(tmp[0], tmp[1], tmp[2]));
		}


		bbCoordList.push_back(new XYZ(bs->ncrd()[0], bs->ncrd()[1], bs->ncrd()[2]));
		bbCoordList.push_back(new XYZ(bs->cacrd()[0], bs->cacrd()[1], bs->cacrd()[2]));
		bbCoordList.push_back(new XYZ(bs->ccrd()[0], bs->ccrd()[1], bs->ccrd()[2]));
		bbCoordList.push_back(new XYZ(bs->ocrd()[0], bs->ocrd()[1], bs->ocrd()[2]));

		string N = triName+"-N";
		string CA = triName+"-CA";
		string C = triName + "-C";
		string O = triName + "-O";
		bbApList.push_back(atLib->getAtomProperty(N));
		bbApList.push_back(atLib->getAtomProperty(CA));
		bbApList.push_back(atLib->getAtomProperty(C));
		bbApList.push_back(atLib->getAtomProperty(O));

		this->hasAromaticRing = rot->hasAromaticRing;
		if(hasAromaticRing){
			this->normalVectorOfRing = cs.local2globalcrd(rot->normalVectorOfRing) - cs.origin_;
		}
	}
	virtual ~Conformer();
};


class ConformerGroup{
	public:
	vector<Conformer*> confList;
	int rotNum;

	ConformerGroup(RotamerGroup* group, BackBoneSite* bs, AtomLib* atLib){
		int n = group->rotNum;
		for(int i=0;i<n;i++){
			confList.push_back(new Conformer(group->rotList.at(i), bs, atLib));
		}
		this->rotNum = n;
	}

	virtual ~ConformerGroup();
};


} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_CONFORMER_H_ */
