/*
 * backbonebuilder.cpp
 *
 *  Created on: 2018年1月9日
 *      Author: hyliu
 */

#include "backbone/backbonebuilder.h"
#include "dstl/randomengine.h"
#include "pdbstatistics/phipsidistr.h"
#include "backbone/closealoop.h"
using namespace NSPproteinrep;
std::vector<BackBoneSite> BackBoneBuilder::buildforwardbackbone(int length,
		const BackBoneSite & nflankingsite, const std::vector<std::pair<int,int>> & helixregions,
		const std::vector<std::pair<int,int>> &strandregions, const std::set<int> & cissites){
	auto rng=NSPdstl::RandomEngine<>::getinstance().realrng(0,1);
	BackBoneSite bsn=nflankingsite;
	std::vector<BackBoneSite> chain(length);
	for(int i=0;i<chain.size();++i){
		bool helix=false;
		bool strand=false;
		for(auto &r:helixregions){
			if(i>=r.first &&i<r.first+r.second) {
				helix=true;break;
			}
		}
		for(auto &r:strandregions){
			if(i>=r.first &&i<r.first+r.second) {
				strand=true;break;
			}
		}
		double phi,psi;
		if(helix) NSPpdbstatistics::PhiPsiDistr::helixdistr().randomphipsi(rng,&phi,&psi);
		else if(strand)NSPpdbstatistics::PhiPsiDistr::stranddistr().randomphipsi(rng,&phi,&psi);
		else NSPpdbstatistics::PhiPsiDistr::mixcoildistr().randomphipsi(rng,&phi,&psi);
		BackBoneSite *bsp=&bsn;
		if(i>0) bsp=&chain[i-1];
		bool cispep=(cissites.find(i) != cissites.end());
		genbackbonesite(bsp,cispep,phi,psi,&chain[i]);
		if(i==0) {
			chain[0].chainid='A';
			chain[0].resseq=0;
		}
		chain[i].resname="ALA";
		if(cispep) chain[i].resname="PRO";
		if(helix) chain[i].sscode='H';
		else chain[i].sscode='C';
	}
	return chain;
}
std::vector<BackBoneSite> BackBoneBuilder::buildbackwardbackbone(int length,
		const BackBoneSite & cflankingsite, const std::vector<std::pair<int,int>> & helixregions,
		const std::vector<std::pair<int,int>> & strandregions, const std::set<int> & cissites){
	auto rng=NSPdstl::RandomEngine<>::getinstance().realrng(0,1);
	BackBoneSite bsc=cflankingsite;
	std::vector<BackBoneSite> chain(length);
	for(int i=chain.size()-1;i>=0;--i){
		bool helix=false;
		bool strand=false;
		for(auto &r:helixregions){
			if(i>=r.first &&i<r.first+r.second) {
				helix=true;break;
			}
		}
		for(auto &r:strandregions){
			if(i>=r.first &&i<r.first+r.second) {
				strand=true;break;
			}
		}
		double phi,psi;
		if(helix) NSPpdbstatistics::PhiPsiDistr::helixdistr().randomphipsi(rng,&phi,&psi);
		else if(strand)NSPpdbstatistics::PhiPsiDistr::stranddistr().randomphipsi(rng,&phi,&psi);
		else NSPpdbstatistics::PhiPsiDistr::mixcoildistr().randomphipsi(rng,&phi,&psi);
		BackBoneSite *bsm=&bsc;
		if(i<chain.size()-1) bsm=&chain[i+1];
		bool cispep=(cissites.find(i) != cissites.end());
		double omiga=180.0;
		if(cispep) omiga=0.0;
		NSPproteinrep::genprevbackbonesite(bsm,omiga,psi,phi,
				&chain[i]);
		if(i==chain.size()-1) {
			chain[i].chainid='A';
			chain[i].resseq=i;
			chain[i].resname="ALA";
		}
		if(helix) chain[i].sscode='H';
		else chain[i].sscode='C';
	}
	return chain;
}

std::vector<BackBoneSite> BackBoneBuilder::buildstrandat(int length,NSPgeometry::XYZ r0,
		NSPgeometry::XYZ direction,bool forward){
	BackBoneSite nsite;
	genbackbonesite(nullptr, false, 0.0,0.0, &nsite);
	std::vector<std::pair<int,int>> strandregions;
	strandregions.push_back(std::make_pair(0,length));
	std::vector<BackBoneSite> chain=buildforwardbackbone(length,nsite, std::vector<std::pair<int,int>>(),strandregions,std::set<int> ());
	movechainto(r0,direction,forward,chain);
	return chain;
}
std::vector<BackBoneSite> BackBoneBuilder::buildhelixat(int length,NSPgeometry::XYZ r0,
		NSPgeometry::XYZ direction,bool forward){
	BackBoneSite nsite;
	genbackbonesite(nullptr, false, 0.0,0.0, &nsite);
	std::vector<std::pair<int,int>> helixregions;
	helixregions.push_back(std::make_pair(0,length));
	std::vector<BackBoneSite> chain=buildforwardbackbone(length,nsite,
			helixregions,std::vector<std::pair<int,int>>(),std::set<int> ());
	movechainto(r0,direction,forward,chain);
	return chain;
}

void BackBoneBuilder::movechainto(NSPgeometry::XYZ r0,NSPgeometry::XYZ direction,
		bool forward, std::vector<BackBoneSite> & chain){
	NSPgeometry::XYZ rhead;
	NSPgeometry::XYZ rheadtotail;
	if (forward){
		rhead=chain[0].cacrd();
		rheadtotail=chain.back().cacrd()-rhead;
	} else {
		rhead=chain.back().cacrd();
		rheadtotail=chain[0].cacrd()-rhead;
	}
	NSPgeometry::XYZ trans=r0-rhead;
	for (auto &bs:chain){
		bs.translate(trans);
	}
	NSPgeometry::XYZ axis=NSPgeometry::cross(rheadtotail,direction);
	if(axis.squarednorm()>0.0){
		double angle=acos(NSPgeometry::dot(rheadtotail,direction)/(sqrt(rheadtotail.squarednorm()*direction.squarednorm())));
		NSPgeometry::Rotation rot(NSPgeometry::QuaternionCrd(axis,angle,1.0),r0);
		for(auto &bs:chain){
			bs.rotate(rot);
		}
	}
}
std::vector<std::shared_ptr<std::vector<BackBoneSite>>> BackBoneBuilder::buildlinkers(int length,
		const BackBoneSite &nflankingsite,
		const BackBoneSite &cflankingsite,
		const std::vector<std::pair<int,int>> & helixregions,
		const std::vector<std::pair<int,int>> & strandregions,
		const std::set<int> & cissites){
		std::vector<BackBoneSite> temploop(length+2);
		temploop[0]=nflankingsite;
		temploop.back()=cflankingsite;
		std::vector<BackBoneSite> middle=buildforwardbackbone(length,nflankingsite, helixregions,strandregions,cissites);
		temploop[1]=middle[0];
		std::vector<BackBoneSite> revmiddle=buildbackwardbackbone(1,cflankingsite, helixregions,strandregions,cissites);
		temploop[length]=revmiddle[0];
		std::vector<double> torsions;
		for(auto & bs:middle) {
			torsions.push_back(bs.phi());
			torsions.push_back(bs.psi());
			torsions.push_back(bs.omiga());
		}
		CloseALoop closer(temploop,1,length,torsions);
		std::set<int> fixedsites;
		for(int i=0;i<length;++i){
			bool fixedi=false;
			for(auto &h:helixregions) if(i>=h.first && i<h.first+h.second) fixedi=true;
			for(auto &s:strandregions) if(i>=s.first && i<s.first+s.second) fixedi=true;
			if(fixedi) fixedsites.insert(i+1);
		}
		std::vector<std::shared_ptr<std::vector<BackBoneSite>>> solutions=closer.getsolutions(fixedsites);
		return solutions;
}


