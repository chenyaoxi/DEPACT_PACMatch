/*
 * chainassembler.cpp
 *
 *  Created on: 2017年1月4日
 *      Author: hyliu
 */

#include "backbone/chainassembler.h"
#include "backbone/backboneloop.h"
#include  "dstl/randomengine.h"
using namespace NSPproteinrep;


bool NSPproteinrep::assembleSSelements(const std::vector<std::vector<BackBoneSite>> & SS,
		const std::vector<unsigned int> &order,
		const std::map<std::pair<unsigned int,unsigned int>,unsigned int> & looplengths){
	ChainAssembler assembler;
	for(auto & ss:SS) assembler.addSSelement(ss);
	assembler.setelementorder(order);
	assembler.setlooplengths(looplengths);
//control parameters
	unsigned int ntries=100u;
	unsigned int maxsteps=100u*(SS.size()-1);
	unsigned int randomseed=37u;
	unsigned int nprint=100;
	std::string output="assembled";
	IdealGeometries::getGlobalInstance("idealgeometries.dat");

	if(!assembler.rebuildallloops(ntries)){
		std::cout <<"first round loop building failed."<<std::endl;
		return false;
	}
	NSPdstl::RandomEngine<> &reng=NSPdstl::RandomEngine<>::getinstance();
	reng.init(randomseed);

	unsigned int nstep=0;
	while(nstep < maxsteps) {
		if(nstep==(nstep/nprint)*nprint) {
			double score;
			std::vector<BackBoneSite> chain=assembler.getassembledchain(&score);
			std::ofstream ofs;
			std::string filename=output+std::to_string(nstep) + ".pdb";
			ofs.open(filename.c_str());
			ofs <<"COMMENT score= " <<score <<std::endl;
			writeSitesToPDB(ofs,chain);
			ofs.close();
		}
		reng.setintrng(0,order.size()-2);
		unsigned int l=reng.randomint();
		assembler.rebuildloop(order[l],order[l+1],ntries);
		++nstep;
		std::cout <<"nstep: " << nstep <<std::endl;
	}
	return true;
}
void ChainAssembler::addSSelement(const Fragment &ele) {
	SSelements_.push_back(ele);
}
void ChainAssembler::setelementorder(const std::vector<unsigned int> & order) {
	assert(order.size() == SSelements_.size());
	elementorder_.clear();
	std::vector<bool> used(order.size(), false);
	for (auto i : order) {
		assert(i < SSelements_.size());
		assert(!used[i]);
		used[i] = true;
		elementorder_.push_back(i);
	}
}
void ChainAssembler::setlooplengths(
		const std::map<std::pair<unsigned int, unsigned int>, unsigned int> & looplengths) {
	loops_.clear();
	for (auto &l : looplengths) {
		assert(l.first.first < SSelements_.size());
		assert(l.first.second < SSelements_.size());
		loops_.insert(std::make_pair(l.first, Loop()));
		loops_.at(l.first).length = l.second;
	}
	for (unsigned int i = 0; i < elementorder_.size() - 1; ++i) {
		assert(
				loops_.find(
						std::make_pair(elementorder_[i], elementorder_[i + 1]))
						!= loops_.end());
	}
}

bool ChainAssembler::rebuildloop(unsigned int headelement,
		unsigned int tailelement, unsigned int ntries) {
	Loop *loop = &(loops_.at(std::make_pair(headelement, tailelement)));
	std::vector<Fragment> context;
	bool contextcomplete { false };
	makeloopcontext(headelement, tailelement, &context, &contextcomplete);
	BackBoneLoop rloop(&context, 0u, 1u, loop->length);
	std::vector<Fragment> newloops(1, Fragment());
	std::vector<double> score;
	if (rloop.getLoops(1u, ntries, newloops.begin(), &score) <= 0)return false;
	if (score[0] >= loop->score)
		return false;
	loop->loop = std::shared_ptr < Fragment > (new Fragment(newloops[0]));
	if (contextcomplete)
		loop->score = score[0];
	return true;
}

bool ChainAssembler::rebuildallloops(unsigned int ntries) {
	bool success=true;
	for(unsigned int i=0; i<elementorder_.size()-1;++i) {
		if(!rebuildloop(elementorder_[i],elementorder_[i+1],ntries)){
			std::cout <<"Not able to (re)build loop between element "
					<<elementorder_[i] << " and " <<elementorder_[i+1] <<std::endl;
			success=false;;
		}
	}
	return success;
}

std::vector<BackBoneSite> ChainAssembler::getassembledchain(double *score) {
	Fragment chain;
	for (unsigned int i = 0; i < elementorder_.size(); ++i) {
		Fragment &element=SSelements_[elementorder_[i]];
		if( i !=elementorder_.size()-1)
			assert(loopbuilt(elementorder_[i],elementorder_[i+1]));
		auto itbegin=element.begin();
		if(i !=0 ) ++itbegin;
		auto itend=element.end();
		if(i !=elementorder_.size()-1) itend= element.begin()+element.size()-1;
		unsigned int oldsize=chain.size();
		unsigned int inc=itend-itbegin;
		chain.resize(oldsize+inc,BackBoneSite());
		std::copy(itbegin,itend, chain.begin()+oldsize);
		if( i!=elementorder_.size()-1) {
			oldsize = chain.size();
			Fragment *nextloop = loops_.at(
						std::make_pair(elementorder_[i],
								elementorder_[i + 1])).loop.get();
			inc = nextloop->size();
			chain.resize(oldsize + inc,
						BackBoneSite());
			std::copy(nextloop->begin(), nextloop->end(),
						chain.begin() + oldsize);
		}
	}
	//todo scorechain
	return std::move(chain);
}
void ChainAssembler::makeloopcontext(unsigned int headelement,
		unsigned tailelement, std::vector<Fragment> *context, bool *complete) {
	bool inhead = true;
	*complete = true;
	context->resize(2, Fragment());
	for (unsigned int i = 0; i < elementorder_.size(); ++i) {
		if (elementorder_[i] == tailelement)
			inhead = false;
		if (inhead) {
			if (elementorder_[i] != headelement) {
				bool previousloopbuilt;
				if (i == 0)
					previousloopbuilt = false;
				else
					previousloopbuilt = loopbuilt(elementorder_[i - 1],
							elementorder_[i]);
				bool nextloopbuilt = loopbuilt(elementorder_[i],
						elementorder_[i] + 1);
				auto itbegin = SSelements_[elementorder_[i]].begin();
				if (previousloopbuilt)
					++itbegin;
				auto itend = SSelements_[elementorder_[i]].end();
				if (nextloopbuilt)
					itend = SSelements_[elementorder_[i]].begin()
							+ SSelements_[elementorder_[i]].size() - 1;
				if (!nextloopbuilt)
					*complete = false;
				unsigned int oldsize = (*context)[0].size();
				unsigned int sizeincrease = itend - itbegin;
				(*context)[0].resize(oldsize + sizeincrease, BackBoneSite());
				std::copy(itbegin, itend, (*context)[0].begin() + oldsize);
				if (nextloopbuilt) {
					oldsize = (*context)[0].size();
					Fragment *nextloop = loops_.at(
							std::make_pair(elementorder_[i],
									elementorder_[i + 1])).loop.get();
					sizeincrease = nextloop->size();
					(*context)[0].resize(oldsize + sizeincrease,
							BackBoneSite());
					std::copy(nextloop->begin(), nextloop->end(),
							(*context)[0].begin() + oldsize);
				}
			} else {
				bool previousloopbuilt;
				if (i == 0)
					previousloopbuilt = false;
				else
					previousloopbuilt = loopbuilt(elementorder_[i - 1],
							elementorder_[i]);
				auto itbegin = SSelements_[elementorder_[i]].begin();
				if (previousloopbuilt)
					++itbegin;
				auto itend = SSelements_[elementorder_[i]].end();
				unsigned int oldsize = (*context)[0].size();
				unsigned int sizeincrease = itend - itbegin;
				(*context)[0].resize(oldsize + sizeincrease, BackBoneSite());
				std::copy(itbegin, itend, (*context)[0].begin() + oldsize);
			}
		} else {
			if (elementorder_[i] != tailelement) {
				bool previousloopbuilt=loopbuilt(elementorder_[i - 1],
							elementorder_[i]);
				bool nextloopbuilt;
				if(i==elementorder_.size()-1) nextloopbuilt=false;
				else nextloopbuilt= loopbuilt(elementorder_[i],
						elementorder_[i] + 1);
				auto itbegin = SSelements_[elementorder_[i]].begin();
				if (previousloopbuilt)
					++itbegin;
				auto itend = SSelements_[elementorder_[i]].end();
				if (nextloopbuilt)
					itend = SSelements_[elementorder_[i]].begin()
							+ SSelements_[elementorder_[i]].size() - 1;
				if (!nextloopbuilt && i!=elementorder_.size()-1)
					*complete = false;
				unsigned int oldsize = (*context)[1].size();
				unsigned int sizeincrease = itend - itbegin;
				(*context)[1].resize(oldsize + sizeincrease, BackBoneSite());
				std::copy(itbegin, itend, (*context)[1].begin() + oldsize);
				if (nextloopbuilt) {
					oldsize = (*context)[1].size();
					Fragment *nextloop = loops_.at(
							std::make_pair(elementorder_[i],
									elementorder_[i + 1])).loop.get();
					sizeincrease = nextloop->size();
					(*context)[1].resize(oldsize + sizeincrease,
							BackBoneSite());
					std::copy(nextloop->begin(), nextloop->end(),
							(*context)[1].begin() + oldsize);
				}
			} else {
				bool nextloopbuilt;
				if (i == elementorder_.size()-1)
					nextloopbuilt = false;
				else{
					nextloopbuilt = loopbuilt(elementorder_[i],
							elementorder_[i+1]);
					if( !nextloopbuilt) *complete=false;
				}
				auto itbegin = SSelements_[elementorder_[i]].begin();
				auto itend = SSelements_[elementorder_[i]].end();
				if (nextloopbuilt)
						itend = SSelements_[elementorder_[i]].begin()
								+ SSelements_[elementorder_[i]].size() - 1;
				unsigned int oldsize = (*context)[1].size();
				unsigned int sizeincrease = itend - itbegin;
				(*context)[1].resize(oldsize + sizeincrease, BackBoneSite());
				std::copy(itbegin, itend, (*context)[1].begin() + oldsize);
			}
		}
	}
}
