/*
 * backbonetrailmove.cpp
 *
 *  Created on: 2017年7月28日
 *      Author: hyliu
 */


#include "backbone/backbonetrialmoves.h"
#include "dstl/randomengine.h"
using namespace NSPproteinrep;
std::pair<int,int> BackBoneTrialMoves::ChangeRange::chooseposi(int lmin, int lmax,
		int pmin,int pmax){
	NSPdstl::RandomEngine<> &rng=NSPdstl::RandomEngine<>::getinstance();
	int pstart=rng.intrng(pmin,pmax)();
	int len=rng.intrng(lmin,lmax)();
	return std::make_pair(pstart,len);
}
std::pair<int,int> BackBoneTrialMoves::ChangeRange::chooseposi(int chainlength) const{
	int pmax=posimax;
	if(pmax<0) pmax=chainlength-1;
	return chooseposi(lengthmin,lengthmax,posimin,pmax);
}
void BackBoneTrialMoves::chooselocation(const BackBoneTrialMoves::ChangeRange *chgrng){
	if(chgrng) setchangerange(*chgrng);
	std::pair<int,int> pl=chgrange_.chooseposi(chainlength_);
	loopstart_=pl.first;
	looplength_=pl.second;
}
std::vector<std::shared_ptr<Loop>> & BackBoneTrialMoves::maketrialloops(double stepsize){

}
