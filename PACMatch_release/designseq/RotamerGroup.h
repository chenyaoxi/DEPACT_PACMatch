/*
 * RotamerGroup.h
 *
 *  Created on: 2017Äê10ÔÂ23ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_ROTAMERGROUP_H_
#define DESIGNSEQ_ROTAMERGROUP_H_

#include <vector>
#include <list>
#include "designseq/Rotamer.h"

namespace NSPdesignseq {

using namespace std;

class RotamerGroup {
public:
	vector<Rotamer*> rotList;
	vector<vector<Rotamer*>*> aaRots;
	int rotNum;
	RotamerGroup();
	void clear();
	void addRotamer(Rotamer* rot);
	void deleteRotamer(int i);

	virtual ~RotamerGroup();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_ROTAMERGROUP_H_ */
