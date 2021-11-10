/*
 * backbonesasa.cpp
 *
 *  Created on: 2016年12月10日
 *      Author: hyliu
 */

#include "backbone/backbonesasa.h"

using namespace NSPproteinrep;
using namespace NSPgeometry;
std::vector<NSPgeometry::SASAParameter> BackBoneSASA::backboneSASAparameters{
	{1.4,1.45,256},{1.4,1.5,256},{1.4,1.45,256},{1.4,1.45,128},{1.4,1.5,256},
	{1.4,1.5,256}};

BackBoneSASA::BackBoneSASA(const BackBoneSite & bs){
	init(bs);
}

void BackBoneSASA::init(const BackBoneSite &bs) {
	atomsasa_.clear();
	atomsasa_.push_back(AtomSASA(backboneSASAparameters[N]));
	atomsasa_.push_back(AtomSASA(backboneSASAparameters[CA]));
	atomsasa_.push_back(AtomSASA(backboneSASAparameters[C]));
	atomsasa_.push_back(AtomSASA(backboneSASAparameters[O]));
	if(bs.resname != "GLY" ) atomsasa_.push_back(AtomSASA(backboneSASAparameters[CB]));
	if(bs.resname == "PRO") atomsasa_.push_back(AtomSASA(backboneSASAparameters[CD_PRO]));
	reset(bs);
}

void BackBoneSASA::reset(const BackBoneSite &bs){
	XYZ ncrd=bs.getcrd(BackBoneSite::NCRD);
	XYZ cacrd=bs.getcrd(BackBoneSite::CACRD);
	XYZ ccrd=bs.getcrd(BackBoneSite::CCRD);
	XYZ ocrd=bs.getcrd(BackBoneSite::OCRD);
	atomsasa_[N].init(make_localframe(ncrd,cacrd,ccrd));
	atomsasa_[CA].init(make_localframe(cacrd,ncrd,ccrd));
	atomsasa_[C].init(make_localframe(ccrd,cacrd,ocrd));
	atomsasa_[O].init(make_localframe(ocrd,ccrd,cacrd));
	if(bs.resname != "GLY") {
		XYZ cbcrd=bs.cbcrd(1.53);
		atomsasa_[CB].init(make_localframe(cbcrd,cacrd,ncrd));
	}
	if(bs.resname == "PRO") {
			XYZ cdpro=bs.cd_procrd();
			atomsasa_[CB].init(make_localframe(cdpro,ncrd,cacrd));
	}
	for(unsigned int i=0;i<atomsasa_.size()-1;++i){
		for (unsigned int j=i+1;j<atomsasa_.size();++j) {
			updateAtomSASA(atomsasa_[i],atomsasa_[j]);
		}
	}
}

void NSPproteinrep::updateSASA(BackBoneSASA &s1, BackBoneSASA &s2) {
	for( auto it1=s1.atomsasa_.begin(); it1 !=s1.atomsasa_.end();++it1) {
		for (auto it2= s2.atomsasa_.begin();it2 !=s2.atomsasa_.end();++it2) {
			updateAtomSASA(*it1, *it2);
		}
	}
}

void NSPproteinrep::segmentSASA(const std::vector<BackBoneSite> &segment, std::vector<BackBoneSASA> &sasa) {
	sasa.clear();
	for(const auto & s:segment) {
		sasa.push_back(BackBoneSASA(s));
	}
	for(unsigned int i=0;i< sasa.size()-1;++i) {
		for(unsigned int j=i+1; j<sasa.size();++j) {
			updateSASA(sasa[i],sasa[j]);
		}
	}
}
void NSPproteinrep::updatesegmentSASA(std::vector<BackBoneSASA> &sasa1, std::vector<BackBoneSASA> &sasa2){
	for(auto it1=sasa1.begin();it1 != sasa1.end();++it1){
		for(auto it2=sasa2.begin();it2 !=sasa2.end();++it2) {
			updateSASA(*it1,*it2);
		}
	}
}
