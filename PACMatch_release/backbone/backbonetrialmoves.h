/*
 * backbonetrialmoves.h
 *
 *  Created on: 2017年7月27日
 *      Author: hyliu
 */

#ifndef BACKBONE_BACKBONETRIALMOVES_H_
#define BACKBONE_BACKBONETRIALMOVES_H_
#include "backbone/closealoop.h"

namespace NSPproteinrep {

class BackBoneMoves{
public:
	int startposi;
	int endposi;
	std::vector<std::shared_ptr<Loop>> movedloops;
	const std::vector<BackBoneSite> *hostchain;
	Loop *movedloop(int idx) {return movedloops[idx].get();}

};
class BackboneMoveSelector{
public:
	int minposi{0};
	int maxposi{-1};
	int minlength{3};
	int maxlength{7};
	void selectmoves(const std::vector<BackBoneSite> &hostchain,BackBoneMoves *results);
private:
	void selectlocation(BackBoneMoves *moves) const;
	void buildmoves(BackBoneMoves *moves) const;
};
}



#endif /* BACKBONE_BACKBONETRIALMOVES_H_ */
