/*
 * backbonealignment.cpp
 *
 *  Created on: 2017年8月20日
 *      Author: hyliu
 */
#include "backbone/backbonealignment.h"
#include "geometry/structalign.h"
#include <map>
using namespace NSPproteinrep;
using namespace NSPgeometry;
static std::vector<int> seedtoindices(int seed,int chainsize){
	auto ij=dualindex(seed,chainsize);
	std::vector<int> res;
	res.push_back(ij.first-1);
	res.push_back(ij.first);
	res.push_back(ij.first+1);
	res.push_back(ij.second-1);
	res.push_back(ij.second);
	res.push_back(ij.second+1);
	return res;
}
NSPdstl::AlignedPositions NSPproteinrep::alignsites(const std::vector<BackBoneSite> &seta, const std::vector<BackBoneSite> &setb,
		double rcut){
	double rcut2=rcut*rcut;
	 auto matcher=[rcut2](const BackBoneSite &s1, const BackBoneSite &s2)->bool{
		double raa2=(s1.cacrd()-s2.cacrd()).squarednorm();
		if( raa2>rcut2) return false;
		double rnn2=(s1.ncrd()-s2.ncrd()).squarednorm();
		double rcc2=(s1.ccrd()-s2.ccrd()).squarednorm();
		double rnc2=(s1.ncrd()-s2.ccrd()).squarednorm();
		double rcn2=(s1.ccrd()-s2.ncrd()).squarednorm();
		return (rnn2>rcc2?rnn2:rcc2) <(rnc2<rcn2?rnc2:rcn2);
	};
	return NSPdstl::SetMatch::alignset(seta,setb,matcher);
}
std::shared_ptr<BackBoneAlignment::Results> BackBoneAlignment::align(
		const std::vector<BackBoneSite> &confa,
		const std::vector<BackBoneSite> &confb, int seedmode){
	BackBoneAlignment ba(seedmode);
	ba.init(confa,confb);
	std::shared_ptr<Results> optres=std::shared_ptr<Results> (new Results);
	ba.optalign(optres.get());
	return optres;
}
void BackBoneAlignment::init(const std::vector<BackBoneSite> &chaina, const std::vector<BackBoneSite> &chainb){
	assert(chaina.size()>5 &&chainb.size()>5);
	confa_=&chaina;
	confb_=&chainb;
	crda_.clear();
	for(auto &s:chaina) crda_.push_back(s.cacrd());
	for(auto &s:chainb) crdb_.push_back(s.cacrd());
	selectseeds();
	treatedseedb_.clear();
	treatedseedb_.resize(seedsa_.size(),std::set<int>());
}
void BackBoneAlignment::selectseeds(){
	if(seedmode_==AUTO_SEEDS || seedmode_==SS_SEEDS){
		std::vector<int> ssseqa;
		std::vector<int> ssseqb;
		std::vector<std::pair<int,int>> aelmt;
		std::vector<std::pair<int,int>> belmt;
		int nssa=sselements(*confa_,&ssseqa, &aelmt);
		int nssb=sselements(*confb_,&ssseqb,&belmt);
		if((seedmode_==SS_SEEDS && nssa>1 && nssb>1) ||
			(nssa>2 &&nssb>2)){
			SS_selectseeds(ssseqa,ssseqb);
		} else seedmode_=SIMPLE_SEEDS;
	}
	if(seedmode_==SIMPLE_SEEDS){
		simple_selectseeds();
	}
}
void BackBoneAlignment::SS_selectseeds(const std::vector<int> & ssseqa, const std::vector<int> &ssseqb){
	for(int i=1;i<crda_.size()-2;++i) {
		if(ssseqa[i]<0) continue;
		for(int j=i+1;j<crda_.size()-1;++j){
			if(ssseqa[j]<0) continue;
			if(ssseqa[i]==ssseqa[j]) continue;
			double dist2=(crda_[i]-crda_[j]).squarednorm();
			if(dist2>64.0) continue;
			seedsa_.push_back(singleindex(i,j,crda_.size()));
		}
	}
	for(int i=1;i<crdb_.size()-2;++i) {
		if(ssseqb[i] < 0) continue;
		for(int j=i+1;j<crdb_.size()-1;++j){
			if(ssseqb[j]<0) continue;
			if(ssseqb[i]==ssseqb[j]) continue;
			double dist2=(crdb_[i]-crdb_[j]).squarednorm();
			if(dist2>64.0) continue;
			seedsb_.push_back(singleindex(i,j,crdb_.size()));
			seedsb_.push_back(singleindex(j,i,crdb_.size()));
		}
	}
}
void BackBoneAlignment::simple_selectseeds(){
	for(int i=1;i<crda_.size()-4;++i) {
		for(int j=i+3;j<crda_.size()-1;++j){
			double dist2=(crda_[i]-crda_[j]).squarednorm();
			if(dist2>64.0) continue;
			seedsa_.push_back(singleindex(i,j,crda_.size()));
		}
	}
	for(int i=1;i<crdb_.size()-4;++i) {
		for(int j=i+3;j<crdb_.size()-1;++j){
			double dist2=(crdb_[i]-crdb_[j]).squarednorm();
			if(dist2>64.0) continue;
			seedsb_.push_back(singleindex(i,j,crdb_.size()));
			seedsb_.push_back(singleindex(j,i,crdb_.size()));
		}
	}
}
void BackBoneAlignment::optalign(Results *optres){
	Results tryresults;
	Results  &refresults=*optres;
	optres->alignedpositions.clear();
	int goodlength=confa_->size()<confb_->size()? confa_->size()*0.6:confb_->size()*0.6;
	if(goodlength<50) goodlength=50;
	int idxa=0;
//	std::cout <<"Total seed to try: " <<seedsa_.size()<<std::endl;
	bool done=false;
	for(int seeda:seedsa_){
		if(done) break;
		std::set<int> & treatedb=treatedseedb_.at(idxa++);
//		std::cout<<"Trying seed " <<idxa <<std::endl;
//		int idxb=0;
		for(int seedb:seedsb_){
//			std::cout<<"idxb: " <<idxb++<<std::endl;
			if(treatedb.find(seedb) !=
					treatedb.end()) continue;
			tryalign(seeda,seedb, &tryresults);
//			treatedb.insert(seedb);
			if(Results::better(refresults,tryresults)) {
				refresults=tryresults;
				if(refresults.alignedpositions.size()>=goodlength) {
					done=true;
					break;
				}
			}
		}
	}
	std::vector<BackBoneSite> confbt=*confb_;
	for(auto &s:confbt){
		s.rotate(refresults.rigidtransform->rotation());
		s.translate(refresults.rigidtransform->translation());
	}
	NSPdstl::AlignedPositions &ap=refresults.alignedpositions;
	NSPgeometry::RigidTransform &rt=*(refresults.rigidtransform);
	int oldalign=refresults.alignedpositions.size();
	int cycle=0;
	while(true)
	{
		ap=alignsites(*confa_,confbt,3.0);
		if(ap.size()<= oldalign) break;
		oldalign=ap.size();
		rt=superpose(crda_,crdb_,ap,&(refresults.rmsd2));
		confbt=*confb_;
		for(auto &s:confbt) {
			s.rotate(refresults.rigidtransform->rotation());
			s.translate(refresults.rigidtransform->translation());
		}
	}
}
void BackBoneAlignment::tryalign(int seeda,int seedb,
		Results *res){
	auto aindices=seedtoindices(seeda,confa_->size());
	auto bindices=seedtoindices(seedb,confb_->size());
	NSPdstl::AlignedPositions &ap=res->alignedpositions;
	ap.clear();
	for(int i=0;i<aindices.size();++i) ap.push_back(std::make_pair(aindices[i],bindices[i]));
	res->rigidtransform=std::shared_ptr<RigidTransform>(new RigidTransform);
	RigidTransform &rt=*(res->rigidtransform);
	rt=superpose(crda_,crdb_,ap,&(res->rmsd2));
	if(res->rmsd2>3.24) return;
	std::vector<BackBoneSite> confbt=*confb_;
	for(auto &s:confbt){
		s.rotate(rt.rotation());
		s.translate(rt.translation());
	}
	int oldalign=ap.size();
//	int cycle=0;
	while(true)
	{
		ap=alignsites(*confa_,confbt,1.8);
		if(ap.size()<= oldalign) break;
		oldalign=ap.size();
		rt=superpose(crda_,crdb_,ap,&(res->rmsd2));
		confbt=*confb_;
		for(auto &s:confbt) {
			s.rotate(rt.rotation());
			s.translate(rt.translation());
		}
//		if(cycle++ >9) {
//			std::cout<<"superpose cycle: "<< cycle <<std::endl;
//		}
	}
	updatetreated(ap);
	return;
}

void BackBoneAlignment::updatetreated(const NSPdstl::AlignedPositions & ap){
	std::map<int,int> a2b;
	for(auto &a:ap)a2b.insert(a);
	int idxa=0;
	for(auto seeda:seedsa_){
		std::set<int> &treated=treatedseedb_.at(idxa++);
		auto aindices=seedtoindices(seeda,crda_.size());
		std::vector<int> bindices;
		bool update=true;
		for(auto ai:aindices){
			auto it=a2b.find(ai);
			if(it == a2b.end()){
				update=false;
				break;
			}
			bindices.push_back(it->second);
		}
		if(!update) continue;
		if(bindices[1]-bindices[0]!=1 ||
			bindices[2]-bindices[1] !=1 ||
			bindices[4]-bindices[3] !=1 ||
			bindices[5]-bindices[4] !=1) continue;

		int seedb=singleindex(bindices[1],bindices[4],crdb_.size());
		treated.insert(seedb);
	}
}
