/*
 * pepconformer.cpp
 *
 *  Created on: 2017年6月29日
 *      Author: hyliu
 */
#include "backbone/pepconformer.h"

using namespace NSPproteinrep;
using namespace NSPgeometry;
LocalFrame & PepConformer::calclocalframe(){
	int ca=4*(length_/2)+1;
	XYZ rca=globalcrd_[ca];
	XYZ rn=globalcrd_[ca-1];
	XYZ rc=globalcrd_[ca+1];
	localframe_=calclocalframe(rca,rn,rc);
	return localframe_;
}
NSPgeometry::LocalFrame PepConformer::calclocalframe
	(const XYZ &rca, const XYZ &rn, const XYZ &rc){
	XYZ rcb=InternaltoXYZ(rca,rc,rn,1.5,109.5*3.14159265/180.0,
			120*3.14159265/180.0);
	return make_localframe(rca,rcb,0.5*(rc+rn));
}
std::vector<NSPgeometry::XYZ> & PepConformer::calclocalcrd(){
	calclocalframe();
	localcrd_.clear();
	for(auto & r:globalcrd_){
		localcrd_.push_back(localframe_.global2localcrd(r));
	}
	return localcrd_;
}
void PepConformer::tobackbonesites(std::vector<BackBoneSite> &sites,
		const std::vector<NSPgeometry::XYZ> & crd){
	for(int i=0; i<length_;++i) {
		BackBoneSite bs;
		bs.resid=i;
		bs.resseq=i;
		bs.resname="ALA";
		bs.chainid='A';
		std::vector<XYZ> stcrd(4);
		std::copy(crd.begin()+4*i,crd.begin()+4*(i+1),stcrd.begin());
		bs.changecrd(stcrd);
		bs.data_[BackBoneSite::PHI]=360;
		bs.data_[BackBoneSite::PSI]=360;
		bs.data_[BackBoneSite::OMIGA]=360;
		bs.sscode=' ';
		sites.push_back(bs);
	}
}
PepConformer NSPproteinrep::make_conformer(int length,const std::vector<BackBoneSite> & sites, int posi){
	PepConformer res;
	res.length()=length;
	std::vector<XYZ> & globalcrd=res.globalcrd();
	int start=posi-length/2;
	assert(start>=0);
	assert(fragstartsite(sites.begin()+start,sites.end(),length,std::string(),false));
	for(int i=start; i<start+length; ++i){
		globalcrd.push_back(sites[i].getcrd(BackBoneSite::NCRD));
		globalcrd.push_back(sites[i].getcrd(BackBoneSite::CACRD));
		globalcrd.push_back(sites[i].getcrd(BackBoneSite::CCRD));
		globalcrd.push_back(sites[i].getcrd(BackBoneSite::OCRD));
	}
	res.calclocalcrd();
	return res;
}
