/*
 * backbonemoves.h
 *
 *  Created on: 2017年8月1日
 *      Author: hyliu
 */

#ifndef BACKBONE_BACKBONEMOVES_H_
#define BACKBONE_BACKBONEMOVES_H_
#include "backbone/backbonesite.h"
#include "loopclosure/loopsolver.h"
#include <memory>
#include <set>

namespace NSPproteinrep {
typedef std::vector<BackBoneSite> Loop;
class BackBoneMoves{
public:
	int startposi{0};
	int endposi{0};
	std::vector<std::shared_ptr<Loop>> movedloops;
	int nloops() const {return movedloops.size();}
	const std::vector<BackBoneSite> *hostchain;
	Loop *movedloop(int idx) {return movedloops[idx].get();}
	void updatechain(std::vector<BackBoneSite> *chain,bool ssaspbtype,int newloopidx);
	void reinit(){energy0=0.0; energy0set=false;movedloops.clear();}
	double energy0{0.0};
	bool energy0set{false};
};
class BackBoneMoveSelector{
public:
	enum MOVEMODE{UNSET,PERTURBLOOP,SEGMENT};
	enum FIXPOSITIONMODE{NONE,SSELEMENTS,ALLNONGAP,USER};
	BackBoneMoveSelector(){;}
	BackBoneMoveSelector(int chainlength,const std::set<int> &fixedpositions=std::set<int>()){
		setfixedpositions(chainlength,fixedpositions);
	}
	void setfixedpositions(int chainlength,const std::set<int>&fixedpositions);
	int minflexposi{3};
	int maxflexposi{9};
	double maxrotate{10.0};
	double maxtranslate{1.5};
	double torsionsteps{1.0};
	int mode{UNSET};
	void selectmoves(const std::vector<BackBoneSite> &hostchain,
			BackBoneMoves *results){
		if(mode==PERTURBLOOP) perturbloopmoves(hostchain,results);
		else if(mode==SEGMENT) segmentmove(hostchain,results);
		else if(mode==UNSET) {
			std::cout <<" Must set a movemode for backbonemoveselector." << std::endl;
			abort();
		}
	}
	void perturbloopmoves(const std::vector<BackBoneSite> &hostchain,
			BackBoneMoves *results);
	void segmentmove(const std::vector<BackBoneSite> &hostchain,BackBoneMoves *moves);

private:
	void selectlocation(BackBoneMoves *moves);
	void selectlocation(BackBoneMoves *moves,std::set<int> *changepositions);
	std::set<int> fixedpositions_;
	std::vector<int> flexiblelist_;
};
class LoopMover {
public:
	LoopMover(){;}
	LoopMover(int sp,int np):startposi_(sp),endposi_(np){;}
	void makemoves(const std::vector<BackBoneSite> &hostchain, int nnewloops, int ntries,
			BackBoneMoves *result)const;
private:
	int startposi_{-1};
	int endposi_{-1};
};
std::vector<LoopMover> getloopmovers(const std::vector<BackBoneSite> &hostchain,
		int minlooplength,int maxlooplength);
BackBoneMoveSelector make_segmentmoveselector(const std::vector<BackBoneSite> &hostchain,
		int fixmode=BackBoneMoveSelector::NONE,const std::set<int> &fixed=std::set<int>());
BackBoneMoveSelector make_perturbloopselector(const std::vector<BackBoneSite> &hostchain,
		int fixmode=BackBoneMoveSelector::NONE,const std::set<int> &fixed=std::set<int>());
}





#endif /* BACKBONE_BACKBONEMOVES_H_ */
