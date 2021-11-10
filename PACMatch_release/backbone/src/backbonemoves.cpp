/*
 * backbonemoves.cpp
 *
 *  Created on: 2017年8月1日
 *      Author: hyliu
 */
#include "backbone/backbonemoves.h"
#include "dstl/randomengine.h"
#include "backbone/closealoop.h"
#include "pdbstatistics/proteinblock.h"
#include "pdbstatistics/phipsidistr.h"

using namespace NSPproteinrep;
using namespace NSPgeometry;
using namespace NSPdstl;
void BackBoneMoves::updatechain(std::vector<BackBoneSite> *chain,bool ssaspbtype,
		int newloopidx) {
	int length = endposi - startposi;
	int chainlength = hostchain->size();
	if (length <= 0)
		length = chainlength + length;
	const Loop *nl = movedloops[newloopidx].get();
	for (int i = 0; i < length; ++i) {
		int posi = (startposi + i) % chainlength;
		chain->at(posi) = nl->at(i);
	}
	if (startposi < endposi) {
		NSPpdbstatistics::ProteinBlock::setpbtypes(*chain, ssaspbtype,startposi - 2,
				endposi + 2);
	} else {
		NSPpdbstatistics::ProteinBlock::setpbtypes(*chain, ssaspbtype,startposi - 2);
		if (endposi > 0)
			NSPpdbstatistics::ProteinBlock::setpbtypes(*chain, ssaspbtype,0, endposi + 2);
	}
}
static std::set<int> genfixed(const std::vector<BackBoneSite> &hostchain, int fixmode,
		const std::set<int> &fixed){
	std::set<int> myfixed;
	if(fixmode==BackBoneMoveSelector::SSELEMENTS){
		int posi=0;
		for(auto &s:hostchain){
			char ss=s.sscode;
			if(ss=='H' || ss=='E' ||ss=='m' || ss=='d') {
				if(!s.isgap) myfixed.insert(posi);
			}
			posi++;
		}
	} else if(fixmode==BackBoneMoveSelector::ALLNONGAP) {
		int posi=0;
		for(auto &s:hostchain){
			if(!s.isgap) myfixed.insert(posi);
			posi++;
		}
	} else if(fixmode==BackBoneMoveSelector::USER){
		myfixed=fixed;
	}
	for(int i=1;i<hostchain.size();++i){
		if(hostchain[i-1].isgap) myfixed.insert(i);
	}
	return myfixed;
}
BackBoneMoveSelector NSPproteinrep::make_perturbloopselector(const std::vector<BackBoneSite> &hostchain,
		int fixmode,const std::set<int> &fixed){
	BackBoneMoveSelector ms;
	ms.mode=BackBoneMoveSelector::PERTURBLOOP;
	ms.minflexposi=3;
	ms.maxflexposi=7;
	ms.setfixedpositions(hostchain.size(),genfixed(hostchain,fixmode,fixed));
	return ms;
}
BackBoneMoveSelector NSPproteinrep::make_segmentmoveselector(const std::vector<BackBoneSite> &hostchain,
		int fixmode,const std::set<int> &fixed){
		BackBoneMoveSelector ms;
		ms.mode=BackBoneMoveSelector::SEGMENT;
		ms.minflexposi=2;
		std::set<int> myfixed=genfixed(hostchain,fixmode,fixed);
		ms.maxflexposi=(hostchain.size()-myfixed.size()+1)/2;
		ms.setfixedpositions(hostchain.size(),myfixed);
		return ms;
}
void BackBoneMoveSelector::setfixedpositions(int chainlength,
		const std::set<int> & fixed) {
	fixedpositions_ = fixed;
	fixedpositions_.erase(chainlength-1);
	for (int i = 0; i < chainlength; ++i) {
		if (fixed.find(i) == fixed.end())
			flexiblelist_.push_back(i);
	}
	if(maxflexposi>flexiblelist_.size()-1) maxflexposi=flexiblelist_.size()-1;
	assert(maxflexposi>=minflexposi);
	assert(flexiblelist_.size()>minflexposi);
}
void BackBoneMoveSelector::selectlocation(BackBoneMoves *moves) {
	RandomEngine<> & rng = RandomEngine<>::getinstance();
	int chainlength = moves->hostchain->size();
	if (flexiblelist_.empty())
		setfixedpositions(chainlength, std::set<int>());
	int posistart = rng.intrng(1, flexiblelist_.size() - 1)();
	int nposi = rng.intrng(minflexposi, maxflexposi)();
	assert(flexiblelist_.size() >= nposi);
	int posiend = (posistart + nposi-1) % flexiblelist_.size();
	moves->startposi = flexiblelist_[posistart];
	moves->endposi = (flexiblelist_[posiend]+1)%chainlength;
}
void BackBoneMoveSelector::selectlocation(BackBoneMoves *moves,std::set<int> *changepositions) {
	RandomEngine<> & rng = RandomEngine<>::getinstance();
	int chainlength = moves->hostchain->size();
	if (flexiblelist_.empty())
		setfixedpositions(chainlength, std::set<int>());
	int posistart = rng.intrng(0, flexiblelist_.size() - 1)();
	int nposi = rng.intrng(minflexposi, maxflexposi)();
	int nintmax=flexiblelist_.size()/(2*nposi);
	if(nintmax<1)nintmax=1;
	if (nintmax>5) nintmax=5;
	int length=1;
	changepositions->clear();
	for(int i=0;i<nposi;++i){
		changepositions->insert(flexiblelist_[(posistart+length-1)%flexiblelist_.size()]);
		if(i!=nposi-1) length += rng.intrng(1,nintmax)();
	}
	assert(length<=flexiblelist_.size()-2);
	moves->startposi = flexiblelist_[posistart];
	moves->endposi = (flexiblelist_[(posistart+length-1)%flexiblelist_.size()]+1)%chainlength;
}
void BackBoneMoveSelector::perturbloopmoves(
		const std::vector<BackBoneSite> &hostchain, BackBoneMoves *results) {
	typename RandomEngine<>::real_generator_type & rng =
			RandomEngine<>::getinstance().realrng(0.0, 1.0);
	results->reinit();
	results->hostchain = &hostchain;
	std::set<int> changepositions;
//	selectlocation(results);
	bool done=false;
	int chainlength= hostchain.size();
	int length;
	while(!done){
		selectlocation(results,&changepositions);
		if(results->startposi==0) continue;
		length = results->endposi - results->startposi;
		if (length <= 0)
			length = chainlength + length;
		bool containgap=false;
		for(int i=0;i<length;++i){
			int posi = (results->startposi + i) % chainlength;
			if(hostchain[posi].isgap) {
				containgap=true;
				break;
			}
		}
		if(containgap) continue;
		done=true;
	}
	std::vector<double> newtorsions;
	for (int i = 0; i < length; ++i) {
		int posi = (results->startposi + i) % chainlength;
		newtorsions.push_back(hostchain[posi].phi());
		newtorsions.push_back(hostchain[posi].psi());
		newtorsions.push_back(hostchain[posi].omiga());
	}
	int tidx = 0;
	for (int i = 0; i < length; ++i) {
		int posi = (results->startposi + i) % chainlength;
//		if (fixedpositions_.find(posi) == fixedpositions_.end()) {
		if (changepositions.find(posi) != changepositions.end()) {
			newtorsions[tidx] += (2.0 * rng() - 1.0) * torsionsteps;
			newtorsions[tidx + 1] += (2.0 * rng() - 1.0) * torsionsteps;
		}
		tidx += 3;
	}
	results->movedloops.clear();
	if (results->endposi < results->startposi) {
		results->movedloops.push_back(std::shared_ptr < Loop > (new Loop));
		Loop *l = results->movedloops.back().get();
		l->resize(length);
		if (results->startposi > 0) {
			BackBoneSite site0 = hostchain[results->startposi - 1];
			BackBoneSite *psite = &site0;
			bool cispep = site0.nextpepcis();
			int tidx = 0;
			for (int i = results->startposi; i < hostchain.size(); ++i) {
				BackBoneSite *csite = &((*l)[i - results->startposi]);
				genbackbonesite(psite, cispep, newtorsions[tidx],
						newtorsions[tidx + 1], csite);
				csite->resname = hostchain[i].resname;
				psite = csite;
				tidx += 2;
				double omiga = newtorsions[tidx];
				cispep = omiga > -90.0 && omiga < 90.0;
				++tidx;
			}
		}
		if (results->endposi > 0) {
			BackBoneSite siten = hostchain[results->endposi];
			BackBoneSite *nsite = &siten;
			int tidx = newtorsions.size() - 3;
			for (int i = results->endposi - 1; i >= 0; --i) {
				BackBoneSite *csite = &((*l)[length - (results->endposi - i)]);
				genprevbackbonesite(nsite, newtorsions[tidx + 2],
						newtorsions[tidx + 1], newtorsions[tidx], csite);
				csite->resname = hostchain[i].resname;
				nsite = csite;
				tidx -= 3;
			}
		}
	} else {
		CloseALoop cl(hostchain, results->startposi, length, newtorsions);
		std::set<int> fixed;
		for(int i=results->startposi;i<results->endposi;++i) {
			if(changepositions.find(i) == changepositions.end()) fixed.insert(i);
		}
		results->movedloops = cl.getsolutions(fixed);
	}
}
void LoopMover::makemoves(const std::vector<BackBoneSite> &hostchain,
		int nnewloop, int ntries, BackBoneMoves *results) const {
	typename RandomEngine<>::real_generator_type & rng =
			RandomEngine<>::getinstance().realrng(0.0, 1.0);
	results->reinit();
	results->hostchain = &hostchain;
	results->startposi = startposi_;
	results->endposi = endposi_;
	int chainlength = hostchain.size();
	assert(
			startposi_ < endposi_ && startposi_ > 0
					&& endposi_ < chainlength);
	int length = results->endposi - results->startposi;
	int nloops = 0;
	int ntry = 0;
	results->movedloops.clear();
	while (nloops < nnewloop && ntry < ntries) {
		std::vector<double> newtorsions;
		for (int i = 0; i < length; ++i) {
			int posi = (results->startposi + i) % chainlength;
			auto &distr = NSPpdbstatistics::PhiPsiDistr::phipsidistr(
					hostchain[posi].resname, hostchain[posi + 1].resname);
			double phi, psi;
/*			distr.randomphipsi(rng, &phi, &psi);
			while(true) {
				phi=rng()*360.0;
				psi=rng()*360.0;
				double rp=distr.distr(phi,psi);
				if(rng()<rp/1.e-4) break;
			}
			newtorsions.push_back(phi);
			newtorsions.push_back(psi);
			newtorsions.push_back(hostchain[posi].omiga());*/
			newtorsions.push_back(hostchain[posi].phi()+2.0*(rng()-1.0)*5.0);
			newtorsions.push_back(hostchain[posi].psi()+2.0*(rng()-1.0)*5.0);
			newtorsions.push_back(hostchain[posi].omiga());
		}
		CloseALoop cl(hostchain, results->startposi, length, newtorsions);
		auto solutions = cl.getsolutions();
		for (auto &s : solutions)
			results->movedloops.push_back(s);
		nloops += solutions.size();
		++ntry;
	}
}
static RigidTransform gentransform(const BackBoneSite &s1, const BackBoneSite &s2,
		double maxrotate, double maxtranslate){
	static auto &reng=NSPdstl::RandomEngine<>::getinstance();
	double angle=maxrotate*(reng.realrng(0.0,1.0)()-0.5);
	if(!s1.isgap && !s2.isgap) {
		XYZ center=s1.cacrd();
		XYZ axis=s2.cacrd()-center;
		return RigidTransform(QuaternionCrd(axis,angle),center,XYZ(0.0,0.0,0.0));
	} else if(s1.isgap &&!s2.isgap) {
		XYZ center=s2.cacrd();
		XYZ axis=s2.ncrd()-s2.cacrd();
		return RigidTransform(QuaternionCrd(axis,angle),center,XYZ(0.0,0.0,0.0));
	} else if(!s1.isgap && s2.isgap){
		XYZ center=s1.cacrd();
			XYZ axis=s1.ccrd()-s1.cacrd();
			return RigidTransform(QuaternionCrd(axis,angle),center,XYZ(0.0,0.0,0.0));
	} else {
		XYZ axis(reng.realrng(0.0,1.0),1.0);
		return RigidTransform(QuaternionCrd(axis,angle),XYZ(0.0,0.0,0.0),
				XYZ(reng.realrng(0.0,1.0),maxtranslate));
	}
}

void BackBoneMoveSelector::segmentmove(const std::vector<BackBoneSite> &hostchain,BackBoneMoves *moves){
	static auto &reng=NSPdstl::RandomEngine<>::getinstance();
	int posistart;
	int posiend;
	bool done=false;
	int sh=(flexiblelist_.size()+1)/2;
	while (!done){
		posistart = reng.intrng(0, flexiblelist_.size()-minflexposi)();
		int nposi = reng.intrng(minflexposi, maxflexposi)();
		if(posistart+nposi>flexiblelist_.size()) nposi=flexiblelist_.size()-posistart;
		if(nposi >sh) nposi=sh;
		posiend=(flexiblelist_[posistart+nposi-1]+1)%hostchain.size();
		posistart=flexiblelist_[posistart];
		if(posistart==posiend) continue;
		bool allgap=true;
		for(int i=posistart;i<=posiend; ++i) {
			if(!hostchain[i].isgap){
				allgap=false;
				break;
			}
		}
		if(allgap) continue;
		done=true;
	}
	moves->reinit();
	moves->hostchain=&hostchain;
	moves->startposi=posistart;
	moves->endposi=posiend;
	BackBoneSite s1=hostchain[posistart];
	BackBoneSite s2=hostchain[posiend-1];
	if(posistart == 0) {
		s1.isgap=true;
	}
	if(posiend == hostchain.size()){
		s2.isgap=true;
	}
	RigidTransform rt=gentransform(s1,s2,maxrotate,maxtranslate);
	moves->movedloops.clear();
	moves->movedloops.push_back(std::shared_ptr<Loop>(new Loop));
	Loop * segment=moves->movedloop(0);
	segment->resize(posiend-posistart);
	std::copy(hostchain.begin()+posistart, hostchain.begin()+posiend, segment->begin());
	if(s1.isgap &&s2.isgap) {
		double count=0.0;
		XYZ center(0.0,0.0,0.0);
		for(auto &s:*segment){
			if(!s.isgap) {
				count +=1.0;
				center =center+s.cacrd();
			}
		}
		center =-(1.0/count)*center;
		for(auto &s:*segment) {
			s.translate(center);
			s.rotate(rt.rotation());
			s.translate(rt.translation()-center);
		}
	} else {
		for(int p=0;p<segment->size();++p){
			auto &s=segment->at(p);
			if(p==0 && !s1.isgap) {
					XYZ ccrd=s.ccrd();
					rt.apply(&ccrd);
					s.data_[BackBoneSite::CCRD]=ccrd.x_;
					s.data_[BackBoneSite::CCRD+1]=ccrd.y_;
					s.data_[BackBoneSite::CCRD+2]=ccrd.z_;
					XYZ ocrd=s.ocrd();
					rt.apply(&ocrd);
					s.data_[BackBoneSite::OCRD]=ocrd.x_;
					s.data_[BackBoneSite::OCRD+1]=ocrd.y_;
					s.data_[BackBoneSite::OCRD+2]=ocrd.z_;
			} else if(p==segment->size()-1 && !s2.isgap){
					XYZ ncrd=s.ncrd();
					rt.apply(&ncrd);
					s.data_[BackBoneSite::NCRD]=ncrd.x_;
					s.data_[BackBoneSite::NCRD+1]=ncrd.y_;
					s.data_[BackBoneSite::NCRD+2]=ncrd.z_;
			} else
				s.rotate(rt.rotation());
		}
	}
	if(!s1.isgap) {
		(*segment)[0].data_[BackBoneSite::PHI]=(*segment)[0].phi(hostchain[posistart-1]);
		(*segment)[0].data_[BackBoneSite::PSI]=(*segment)[0].psi((*segment)[1]);
	}
	if(!s2.isgap){
		segment->back().data_[BackBoneSite::PHI]=segment->back().phi((*segment)[segment->size()-2]);
		segment->back().data_[BackBoneSite::PSI]=segment->back().psi(hostchain[posiend]);
	}
}
std::vector<LoopMover> NSPproteinrep::getloopmovers(const std::vector<BackBoneSite> &chain,
		int minlooplength,int maxlooplength){
	std::vector<LoopMover>movers;
	for(int i=1;i<chain.size()-minlooplength-1;++i){
		if(chain[i-1].sscode != 'H' && chain[i-1].sscode !='E')continue;
		if(chain[i].sscode !='H' && chain[i].sscode != 'E') continue;
		if(chain[i+1].sscodechar() !='C') continue;
		int startposi=i;
		int endposi=i+1;
		while(endposi <chain.size()-2 && chain[endposi].sscodechar()=='C') ++endposi;
		if(chain[endposi].sscodechar() =='C') continue;
		++endposi;
		if(endposi-startposi <minlooplength || endposi-startposi>maxlooplength) continue;
		movers.push_back(LoopMover(startposi,endposi));
	}
	return movers;
}

