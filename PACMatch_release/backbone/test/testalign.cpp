/*
 * testalign.cpp
 *
 *  Created on: 2017年8月20日
 *      Author: hyliu
 */
#include "backbone/backbonesite.h"
#include "geometry/rigidtransform.h"
#include "backbone/backbonealignment.h"
#include <algorithm>

using namespace NSPproteinrep;

int main(int argc, char **argv) {
	std::vector<BackBoneSite> sitesa;
	readbackbonesites(std::string(argv[1]),sitesa);
//	std::vector<BackBoneSite> sitesb(sitesa.size()-11);
//	std::copy(sitesa.begin()+61,sitesa.begin()+81,sitesb.begin());
//	std::copy(sitesa.begin()+5,sitesa.begin()+55,sitesb.begin()+20);
//	std::copy(sitesa.begin()+81,sitesa.end(),sitesb.begin()+70);
	std::vector<BackBoneSite> sitesb;
	readbackbonesites(std::string(argv[2]),sitesb);
	auto rt=NSPgeometry::randomrigidtransform(5,1.0);
	for(auto &s:sitesb){
		s.rotate(rt.rotation());
		s.translate(rt.translation());
	}
/*	auto aligned=alignaspointsets(sitesa,sitesb);
	for(auto pp:aligned){
		std::cout<<pp.first<<":" <<sitesa[pp.first].cacrd().toString()
				<<"\t"<<pp.second<<sitesb[pp.second].cacrd().toString()
				<<"\t"<<(sitesa[pp.first].cacrd()-sitesb[pp.second].cacrd()).squarednorm()<<std::endl;
	}
	*/
	auto align=BackBoneAlignment::align(sitesa,sitesb);
	NSPgeometry::RigidTransform &rtf=*(align->rigidtransform);
	for(auto pp:align->alignedpositions){
		NSPgeometry::XYZ nca=rtf.applytoCopy(sitesb[pp.second].cacrd());
		std::cout<<pp.first<<":" <<sitesa[pp.first].cacrd().toString()
				<<"\t"<<pp.second<<nca.toString()
				<<"\t"<<(sitesa[pp.first].cacrd()-nca).squarednorm()<<std::endl;
	}
}



