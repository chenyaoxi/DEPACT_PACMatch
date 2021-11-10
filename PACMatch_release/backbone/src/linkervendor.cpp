/*
 * linkervendor.cpp
 *
 *  Created on: 2016年11月24日
 *      Author: hyliu
 */

#include <backbone/linkervendor.h>
#include <cassert>
using namespace NSPproteinrep;
unsigned int LinkerVendor::MAXSEGMENTS {100};
unsigned int LinkerVendor::MAXCOPIES{500};
unsigned int LinkerVendor::SegmentKey(const std::vector<BackBoneSite> *sg){
	auto it=segmentkeys_.find(sg);
	if(it != segmentkeys_.end()) return it->second;
	unsigned int id=segmentheads_.size();
	assert(id <MAXSEGMENTS);
	segmentkeys_.insert(std::make_pair(sg,id));
	segmentheads_.push_back(std::vector<BackBoneSite>());
	std::vector<BackBoneSite> & head=segmentheads_.back();
	segmenttails_.push_back(std::vector<BackBoneSite>());
	std::vector<BackBoneSite> & tail=segmenttails_.back();
	unsigned int l=sg->size();
	assert(l>1);
	head.push_back((*sg)[0]);
//	head.push_back((*sg)[1]);
	tail.push_back((*sg)[l-1]);
	return id;
}

unsigned int LinkerVendor::NewSegmentKey(const std::vector<BackBoneSite> *sg){
	auto it=segmentkeys_.find(sg);
	if(it == segmentkeys_.end()) return SegmentKey(sg);
	eraseLinkers(it->second);
	std::vector<BackBoneSite> & head=segmentheads_[it->second];
	std::vector<BackBoneSite> & tail=segmentheads_[it->second];
	unsigned int l=sg->size();
	assert(l>1);
	head.push_back((*sg)[0]);
	head.push_back((*sg)[1]);
	tail.push_back((*sg)[l-1]);
	return it->second;
}
unsigned int LinkerVendor::LinkerKey(unsigned int sg1key, unsigned int sg2key, unsigned int length){
	assert(sg1key < segmentkeys_.size() && sg2key <segmentkeys_.size());
	return genlinkerkey(sg1key,sg2key,length);
}
bool LinkerVendor::makeLinkers(unsigned int linkerkey, std::vector<std::vector<BackBoneSite>> *linker,
		unsigned int ncandidates,
		unsigned int ncopies,
		bool regenerate){
	auto it=linkers_.find(linkerkey);
	if( it == linkers_.end()) {
		insertLinker(linkerkey);
		linkerstorage_.insert(std::make_pair(linkerkey,std::vector<std::vector<BackBoneSite>>()));
	}
	unsigned int newcopies=ncopies;
	linker->clear();
	linker->resize(ncopies,std::vector<BackBoneSite>());
	unsigned int oldcopies=0;
	if(!regenerate) {
		oldcopies=ncopies < linkerstorage_[linkerkey].size()? ncopies:linkerstorage_[linkerkey].size();
		newcopies=ncopies-oldcopies;
		for(unsigned int i=0; i<oldcopies;++i)
			(*linker)[i]=linkerstorage_[linkerkey][i];
	}
	unsigned int ngenerated=0;
	unsigned ntry=0;
	while (newcopies > 0 && ntry <10){
		ngenerated = ngenerated +
				linkers_[linkerkey].getLinkers(newcopies,ncandidates,linker->begin()+oldcopies+ngenerated);
		newcopies -= ngenerated;
		++ntry;
	}
	if(newcopies > 0) {
		return false;
	}
	if(linkerstorage_[linkerkey].size() <MAXCOPIES)
	for(int i=0; i<ngenerated; ++i) {
		if(i+oldcopies >= MAXCOPIES) break;
		linkerstorage_[linkerkey].push_back((*linker)[oldcopies+i]);
	}
	return true;
}

bool LinkerVendor::makeLoops(unsigned int linkerkey, std::vector<std::vector<BackBoneSite>> *loop,
		unsigned int ncandidates,
		unsigned int ncopies,
		bool regenerate){
	auto it=loops_.find(linkerkey);
	if( it == loops_.end()) {
		insertLoop(linkerkey);
		loopstorage_.insert(std::make_pair(linkerkey,std::vector<std::vector<BackBoneSite>>()));
	}
	unsigned int newcopies=ncopies;
	loop->clear();
	loop->resize(ncopies,std::vector<BackBoneSite>());
	unsigned int oldcopies=0;
	if(!regenerate) {
		oldcopies=ncopies < loopstorage_[linkerkey].size()? ncopies:loopstorage_[linkerkey].size();
		newcopies=ncopies-oldcopies;
		for(unsigned int i=0; i<oldcopies;++i)
			(*loop)[i]=loopstorage_[linkerkey][i];
	}
	unsigned int ngenerated=0;
	unsigned ntry=0;
	while (newcopies > 0 && ntry <10){
		ngenerated = ngenerated +
				loops_[linkerkey].getLoops(newcopies,ncandidates,loop->begin()+oldcopies+ngenerated);
		newcopies -= ngenerated;
		++ntry;
	}
	if(newcopies > 0) {
		return false;
	}
	if(loopstorage_[linkerkey].size() <MAXCOPIES)
	for(int i=0; i<ngenerated; ++i) {
		if(i+oldcopies >= MAXCOPIES) break;
		loopstorage_[linkerkey].push_back((*loop)[oldcopies+i]);
	}
	return true;
}

void LinkerVendor::insertLinker(unsigned int linkerkey){
	unsigned int sg1key=sgkey1(linkerkey);
	unsigned int sg2key=sgkey2(linkerkey);
	assert (sg1key <segmenttails_.size() && sg2key <segmentheads_.size());
	linkers_.insert(std::make_pair(linkerkey,BackBoneLinker(linkerlength(linkerkey),
			segmenttails_[sg1key],segmentheads_[sg2key])));
}
void LinkerVendor::insertLoop(unsigned int linkerkey){
	unsigned int sg1key=sgkey1(linkerkey);
	unsigned int sg2key=sgkey2(linkerkey);
	assert (sg1key <segmenttails_.size() && sg2key <segmentheads_.size());
	loops_.insert(std::make_pair(linkerkey,BackBoneLoop(segmenttails_[sg1key][0],
			segmentheads_[sg2key][0],linkerlength(linkerkey))));
}
void LinkerVendor::eraseLinkers(unsigned int sgkey) {
	std::vector<unsigned int> toerase;
	for(auto & l:linkers_){
		if( sgkey1(l.first) == sgkey || sgkey2(l.first) == sgkey) {
			toerase.push_back(l.first);
		}
	}
	for(auto it=toerase.begin(); it != toerase.end(); ++it){
		linkers_.erase(*it);
		linkerstorage_.erase(*it);
	}
}

