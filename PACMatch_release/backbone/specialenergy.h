/*
 * specialenergy.h
 *
 *  Created on: 2017年8月18日
 *      Author: hyliu
 */

#ifndef BACKBONE_SPECIALENERGY_H_
#define BACKBONE_SPECIALENERGY_H_
#include "backbone/backbonemoves.h"
#include "backbone/chainpack.h"
#include <memory>
namespace NSPproteinrep{
class ChainAggregateEnergy;
class  SpecialEnergy {
public:
	double totalenergy(const ChainPack & cpck,EnergyComponents *ecomp);
	double deltaE(const ChainPack &cpck,ChainPackMoves *moves);
	friend SpecialEnergy makespecialenergy();
protected:
	std::shared_ptr<ChainAggregateEnergy> chainaggregateenergy_;
};

SpecialEnergy makespecialenergy();

}

#endif /* BACKBONE_SPECIALENERGY_H_ */
