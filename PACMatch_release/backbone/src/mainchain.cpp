/*
 * mainchain.cpp
 *
 *  Created on: 2017年4月22日
 *      Author: hyliu
 */

#include "backbone/mainchain.h"
#include "backbone/closingloop.h"
#include "dstl/randomengine.h"


#include <map>
#include <algorithm>
using namespace NSPproteinrep;

void MainChain::write(std::ostream &ofs) {
	ofs<<size()<<std::endl;
	for(auto &site:(*this)){
		ofs <<site.toString();
	}
	if(!rigidsegments_.segments.empty()) {
		ofs <<"NUMBER_RIGID_SEGMENTS "<<rigidsegments_.segments.size() <<std::endl;
		for(auto &s:rigidsegments_.segments) {
			ofs<<"SEGMENT_POSI_LENGTH " <<s.first <<" "<<s.second <<std::endl;
		}
	}
}
bool MainChain::read(std::istream &is) {
	char buffer[120];
	is.getline(buffer,120);
	std::istringstream iss;
	iss.str(std::string(buffer));
	int length;
	iss >> length;
	readbackbonesites(is, length,*this);
	if(size() != length) return false;
	if(is.getline(buffer,120)){
		std::string id;
		iss.str(std::string(buffer));
		iss.seekg(0);
		iss>>id;
		if(id =="NUMBER_RIGID_SEGMENTS") {
			int nseg;
			iss>>nseg;
			for(int i=0; i<nseg;++i){
				is.getline(buffer,120);
				iss.str(std::string(buffer));
				iss.seekg(0);
				int posi,slength;
				iss >>id;
				if(id !="SEGMENT_POSI_LENGTH") return false;
				iss>>posi >>length;
				setrigidbetween(posi, posi+length);
			}
		}
	}
	return true;
}
void MainChain::generaterandom(int length) {
	resize(length,BackBoneSite());
	BackBoneSite *psite;
	psite=nullptr;
	begin()->pdbid="RAND";
	begin()->chainid='A';
	for(int i=0; i<length;++i){
		this->at(i).resname=chooseresiduetype();
	}
	for(int i=0;i<length;++i) {
		BackBoneSite *newsite=&(this->at(i));
		std::string next="ALA";
		if(i<length-1) next=this->at(i+1).resname;
		bool cispep=false;
		double phi;
		double psi;
		choosephipsi(newsite->resname,next,&phi,&psi,&cispep);
		genbackbonesite(psite,cispep,phi,psi,newsite);
		psite=newsite;
	}
}

void MainChain::copysegment(int startposi, int length, MainChain *out, int outposi) const{
	std::copy(begin()+startposi,begin()+startposi+length,out->begin()+outposi);
}
void MainChain::resetresseq(){
	int seq=1;
	for (auto &r:*this) {
		r.chainid='A';
		r.resid=seq;
		r.resseq=seq++;
	}
}
void MainChain::insertgap(int posi, int length) {
	assert(posi <=size());
	assert (!localrigid(posi-1) && !localrigid(posi)); //do not insert in rigid region
	BackBoneSite temp;
	temp.pdbid="NONE";
	temp.chainid='A';
	temp.resid=1;
	temp.resseq=1;
	if(size()>0) {
		temp.pdbid=this->at(posi).pdbid;
	}
	temp.resname="MISSING";
	insert(begin()+posi,length,temp);
	int rid=temp.resid;
	int rseq=temp.resseq;
	for(auto it=begin()+posi+1; it != begin()+posi+length; ++it){
		it->resid=rid++;
		it->resseq=rseq++;
	}
	for(auto it=begin()+posi+length; it !=end(); ++it) {
		it->resid +=length;
		it->resseq +=length;
	}
	for(auto & s:rigidsegments_.segments) {
		if(s.first>=posi) s.first += length;
	}
}
ChangeRegionSelector::ChangeRegionSelector(const MainChain *mc):mc_(mc){
	for(int p=0;p<mc_->size();++p){
		if(mc->localrigid(p)) continue;
		flexiblesites_.push_back(p);
	}
	if(flexiblesites_[0]==0 && flexiblesites_[1]>1)
		flexiblesites_.erase(flexiblesites_.begin());
	if(flexiblesites_.back()== mc->size()-1 &&
			flexiblesites_[flexiblesites_.size()-2]!= mc->size()-2)
		flexiblesites_.pop_back();
}
bool ChangeRegionSelector::select4seg(int *startposi, int *endposi) const {
	NSPdstl::RandomEngine<>& rneg=NSPdstl::RandomEngine<>::getinstance();
	bool success=false;
	if(flexiblesites_.size()<3) {
		return false;
	} else {
		int ntry=0;
		while (ntry<100) {
			++ntry;
			int s=rneg.intrng(0,flexiblesites_.size()-3)();
			*startposi=flexiblesites_[s];
			int lencan=1;
			for(int p=*startposi+1; p<*startposi+4;++p){
				if(p>=mc_->size()) continue;
				if(mc_->localrigid(p)) break;
				++lencan;
			}
			if(lencan < 3) continue;
			if(lencan==3) {
				if(*startposi >0)
					if(!mc_->localrigid(*startposi-1)) continue;
			}
			*endposi=*startposi+lencan;
			success=true;
			break;
		}
	}
	return success;
}
bool ChangeRegionSelector::selectinoneloop(int *startposi, int *endposi) const {
	NSPdstl::RandomEngine<>& rneg=NSPdstl::RandomEngine<>::getinstance();
	bool success=false;
	if(flexiblesites_.size()<3) {
		return false;
	} else {
		int ntry=0;
		while (ntry<100) {
			++ntry;
			int s=rneg.intrng(0,flexiblesites_.size()-3)();
			*startposi=flexiblesites_[s];
			int lencan=1;
			for(int p=*startposi+1; p<mc_->size();++p){
				if(mc_->localrigid(p)) break;
				++lencan;
			}
			if(lencan < 3) continue;
			int w=rneg.intrng(3,lencan)();
			*endposi=*startposi+w;
			success=true;
			break;
		}
	}
	return success;
}
void NSPproteinrep::mainchainfromsegments(const std::vector<std::vector<BackBoneSite>> &segments,
		const std::vector<int> & order,const std::vector<int> &looplengths,MainChain *mc) {
	assert(order.size()==segments.size());
	assert(looplengths.size() == segments.size()-1);
	mc->clear();
	int totallength=0;
	for(int s=0; s<segments.size();++s) totallength+=segments[s].size();
	mc->resize(totallength,BackBoneSite());
	int posi0=0;
	for(int s:order) {
		std::copy(segments[s].begin(),segments[s].end(),mc->begin()+posi0);
		posi0 +=segments[s].size();
	}
	mc->resetresseq();
	int posig=0;
	for(int l=0; l<looplengths.size();++l) {
		posig +=segments[order[l]].size();
		mc->insertgap(posig,looplengths[l]);
		posig += looplengths[l];
	}
	int segposi=0;
	int l=0;
	for(int s:order) {
		mc->setrigidbetween(segposi,segposi+segments[s].size());
		if(l<looplengths.size())segposi += segments[s].size()+looplengths[l++];
	}
}
