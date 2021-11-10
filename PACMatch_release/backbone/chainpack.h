/*
 * chainpack.h
 *
 *  Created on: 2017年8月17日
 *      Author: hyliu
 */

#ifndef BACKBONE_CHAINPACK_H_
#define BACKBONE_CHAINPACK_H_
#include "backbone/backbonesite.h"
#include "backbone/backboneenergy.h"
#include "dstl/symmatrix1d.h"
namespace NSPproteinrep{

class ChainPack:public std::vector<std::vector<BackBoneSite>> {
public:
	ChainPack(){;}
	ChainPack(const std::vector<std::vector<BackBoneSite>> & chains):
		std::vector<std::vector<BackBoneSite>>(chains){;}
	typedef NSPdstl::SymMatrix1D<double> EMatrix;
	ChainEnergyControl &econtrol(int i) const {return econtrols_[i];}
	double & energy(int c1,int c2) const {return ematrix_(c1,c2);}
	double & energy(int idx) const {return ematrix_.at(idx);}
	std::set<int> &fixed() {return fixed_;}
	const std::set<int> &fixed() const {return fixed_;}
	void setecontrol(const ChainEnergyControl &ectl, int idx) {
		assert(idx<this->size());
		econtrols_[idx]=ectl;}
	void initematrix() const {ematrix_.init(this->size());}
	void init(const std::string &eneparaset=std::string()){
		for(auto &chain:*this){
			econtrols_.push_back(prepareenergycontrol(&chain,eneparaset));
		}
		initematrix();
	}
private:
	mutable std::vector<ChainEnergyControl> econtrols_;
	mutable EMatrix ematrix_;
	std::set<int> fixed_;
};
class ChainPackMoves{
public:
	std::set<int> moveset;
	std::vector<std::vector<std::shared_ptr<std::vector<BackBoneSite>>>> movedchains;
	void updatechainpack(ChainPack *cpck,int mvidx=0);
	double & energy(int c1,int c2,int mvidx=0) {return ematrices[mvidx](c1,c2);}
	std::vector<BackBoneSite> &movedchain(int idx,int mvidx=0)
			{ return *(movedchains[idx][mvidx]);}
	std::vector<ChainPack::EMatrix> ematrices;
	std::vector<bool> localmove;
	int nmoves(){return movedchains[*moveset.begin()].size();}
};
class CpckMoveSelector {
public:
	std::set<int> selectmoveset(const std::vector<int> & completeset,int nmove);
	void selectmoveset(const ChainPack & cpck,ChainPackMoves *moves);
	void makemoves(const ChainPack &cpck,ChainPackMoves *moves,bool newmoveset=true);
	double singlechainmoveprob{0.8};
	double maxrotate{5.0};
	double maxtranslate{0.5};
};
}



#endif /* BACKBONE_CHAINPACK_H_ */
