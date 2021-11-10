/*
 * chainpack.cpp
 *
 *  Created on: 2017年8月17日
 *      Author: hyliu
 */

#include "backbone/chainpack.h"
#include "dstl/randomengine.h"
using namespace NSPproteinrep;
void ChainPackMoves::updatechainpack(ChainPack *cpck,int mvidx){
	for(auto c:moveset){
		auto & chain=cpck->at(c);
		std::vector<BackBoneSite> & nchain=*(movedchains[c][mvidx]);
		std::copy(nchain.begin(),nchain.end(),chain.begin());
		for(int i=0;i<cpck->size();++i) {
			if(moveset.find(i) != moveset.end()) continue;
			cpck->energy(c,i)=ematrices[mvidx](c,i);
		}
	}
}
std::set<int> CpckMoveSelector::selectmoveset(const std::vector<int> &completeset,int maxset){
	assert(maxset <=completeset.size());
	auto & reng=NSPdstl::RandomEngine<>::getinstance();
	std::set<int> moveset;
	int nmove=1;
	if(reng.realrng(0.0,1.0)()>singlechainmoveprob &&maxset>1)
		nmove=reng.intrng(2,maxset)();
	while(moveset.size()<nmove) {
		moveset.insert(completeset[reng.intrng(0,completeset.size()-1)()]);
	}
	return moveset;
}
void CpckMoveSelector::selectmoveset(const ChainPack & cpck,ChainPackMoves *moves){
	std::vector<int> allowed;
	for(int i=0;i<cpck.size();++i) {
		if(cpck.fixed().find(i) != cpck.fixed().end()) continue;
		allowed.push_back(i);
	}
	int maxset=cpck.size()/2;
	if(maxset>allowed.size()) maxset=allowed.size();
	moves->moveset=selectmoveset(allowed,maxset);
}
void CpckMoveSelector::makemoves(const ChainPack &cpck,ChainPackMoves *moves,bool newmoveset){
	if(newmoveset) selectmoveset(cpck,moves);
	std::set<int> & moveset=moves->moveset;
	auto & movedchains=moves->movedchains;
	movedchains.clear();
	movedchains.resize(cpck.size());
	moves->localmove.clear();
	moves->localmove.resize(cpck.size(),false);
	moves->ematrices.clear();
	NSPgeometry::XYZ center(0.0,0.0,0.0);
	double nsitesmoved{0.0};
	NSPgeometry::RigidTransform rt=NSPgeometry::randomrigidtransform(maxrotate,maxtranslate);
	for(auto c:moveset){
		auto & chain=cpck.at(c);
		nsitesmoved +=chain.size();
		for(auto &s:chain) {center = center+s.cacrd();}
	}
	center=(1.0/nsitesmoved)*center;
	for(auto c:moveset){
		auto & chain=cpck.at(c);
		movedchains[c].push_back(std::shared_ptr<std::vector<BackBoneSite>>(new std::vector<BackBoneSite>));
		std::vector<BackBoneSite> *nchain=movedchains[c].back().get();
		nchain->resize(chain.size());
		std::copy(chain.begin(),chain.end(),nchain->begin());
		for(auto &s:*nchain){
			s.translate(-center);
			s.rotate(rt.rotation());
			s.translate(rt.translation()+center);
		}
	}
	moves->ematrices.push_back(ChainPack::EMatrix());
}

