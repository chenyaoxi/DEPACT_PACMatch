/*
 * coresegments.cpp
 *
 *  Created on: 2017年5月11日
 *      Author: hyliu
 */

#include "backbone/coresegments.h"
#include "backbone/energyfunctions.h"
#include "backbone/closingloop.h"
//#include "dstl/randomengine.h"
#include "dstl/sortindex.h"
#include "dstl/topn.h"
#include <sstream>
#include <fstream>
using namespace NSPproteinrep;
void checklooptorsions(MainChain &mc,int loopstart, int looplength,
		const std::vector<BackBoneSite> &newloop){
	std::vector<BackBoneSite> loop;
	loop.push_back(mc.at(loopstart-1));
	for(auto &bs:newloop) loop.push_back(bs);
	loop.push_back(mc.at(loopstart+looplength));
	for( int i=1; i<loop.size()-1;++i) {
		BackBoneSite &bs=loop.at(i);
		BackBoneSite &ps=loop.at(i-1);
		BackBoneSite &ns=loop.at(i+1);
		double phi=bs.phi();
		double psi=bs.psi();
		double omiga=ps.omiga();
		double rphi=bs.phi(ps);
		double rpsi=bs.psi(ns);
		double romiga=ps.omiga(bs);
//		std::cout <<i<<"\t"<<phi<<"\t"<<psi <<"\t"<<rphi<<"\t"<<rpsi<<std::endl;
		double diff_p=phi-rphi;
		double diff_psi=psi-rpsi;
		double diff_o=omiga-romiga;
		auto shift=[](double t){if(t>180.0) t-=360.0; if(t<-180.0) t+=360.0;return t;};
		auto large=[](double d,double d2=1.e-5){return d>d2 || d<-d2;};
		if(large(shift(diff_p)) || large(shift(diff_psi))||large(shift(diff_o),10.0)){
			std::cout <<"Inconsistent site "<<diff_p <<"\t"<<diff_psi <<"\t"<<diff_o<<std::endl;
			exit(1);
		}
	}
}
std::string CoreSegments::linkedchainname(const LinkedChain &lc) const {
	std::string nm;
	for(auto s:lc.segments) nm+=std::to_string(s);
	for (int sl = 0; sl < lc.loopidx_.size(); sl++) {
		int ss1 = lc.segments[sl];
		int ss2 = lc.segments[sl + 1];
		auto cl=closedloops_.at(std::make_pair(ss1,ss2))[lc.loopidx_[sl]];
		int len=cl.loop->size()-2-cl.nextension-cl.cextension;
		nm +="_"+std::to_string(len);
	}
	return nm;
}
int CoreSegments::linkedchainlength(const LinkedChain & lc) const{
	int len=0;
	for (int sl = 0; sl < lc.loopidx_.size(); sl++) {
		int ss1 = lc.segments[sl];
		int ss2 = lc.segments[sl + 1];
		auto cl=closedloops_.at(std::make_pair(ss1,ss2))[lc.loopidx_[sl]];
		len += at(ss1).size()-cl.nextension+cl.loop->size()-cl.cextension-2;
		if(sl==lc.loopidx_.size()-1) len += at(ss2).size();
	}
	return len;
}

int CoreSegments::linkedchainlooplength(const LinkedChain &lc) const{
	int len=0;
	for (int sl = 0; sl < lc.loopidx_.size(); sl++) {
		int ss1 = lc.segments[sl];
		int ss2 = lc.segments[sl + 1];
		auto cl=closedloops_.at(std::make_pair(ss1,ss2))[lc.loopidx_[sl]];
		len += cl.loop->size();
	}
	return len;
}

double CoreSegments::linkedchainaverageene(const LinkedChain &lc) const{
	return lc.energy/(double) linkedchainlooplength(lc);
}
void CoreSegments::buildlinkedmainchain(const LinkedChain &lc, MainChain *mc) const{
	std::vector<int> gaps;
	for (int sl = 0; sl < lc.loopidx_.size(); sl++) {
		int ss1 = lc.segments[sl];
		int ss2 = lc.segments[sl + 1];
		auto cl=closedloops_.at(std::make_pair(ss1,ss2))[lc.loopidx_[sl]];
		gaps.push_back(cl.loop->size()-cl.nextension-cl.cextension-2);
	}
	buildmainchain(lc.segments,gaps,mc);
	int gposi=0;
	for (int sl = 0; sl < lc.loopidx_.size(); sl++) {
		int ss1 = lc.segments[sl];
		int ss2 = lc.segments[sl + 1];
		auto cl=closedloops_.at(std::make_pair(ss1,ss2))[lc.loopidx_[sl]];
		const std::vector<BackBoneSite> *newloop=cl.loop.get();
		int loopstart=gposi+at(ss1).size()-cl.nextension-1;
		std::copy(newloop->begin(),newloop->end(),mc->begin()+loopstart);
		gposi += at(ss1).size()+gaps[sl];
	}
	mc->resetresseq();
}
bool CoreSegments::read(std::istream &is) {
	int nsegments;
	std::istringstream iss;
	std::vector<int> segmentsizes;
	char buffer[120];
	is.getline(buffer, 120);
	iss.str(std::string(buffer));
	iss >> nsegments;
	segmentsizes.resize(nsegments, 0);
	std::vector<std::pair<int,int>> flexibles(nsegments,{0,0});
	for (int i = 0; i < nsegments; ++i) {
		is.getline(buffer, 120);
		iss.str(std::string(buffer));
		iss.seekg(0);
		iss >> segmentsizes[i] >>flexibles[i].first >>flexibles[i].second;
	}
	bool success = false;
	clear();
	n_flexible_.clear();
	c_flexible_.clear();
	for (int i = 0; i < nsegments; ++i) {
		push_back(std::vector<BackBoneSite>());
		success = readbackbonesites(is, segmentsizes[i], back());
		if (!success)
			return success;
		n_flexible_.push_back(flexibles[i].first);
		c_flexible_.push_back(flexibles[i].second);
	}
	return success;
}
void CoreSegments::buildmainchain(const std::vector<int> & order,
		const std::vector<int> & gaps, MainChain *mc) const{
	mc->clear();
	int totallength=0;
	assert(order.size()<=size());
	for(int i=0; i<order.size();++i) totallength+=at(order[i]).size();
	mc->resize(totallength,BackBoneSite());
	int posi0=0;
	for(int s:order) {
		std::copy(at(s).begin(),at(s).end(),mc->begin()+posi0);
		posi0 +=at(s).size();
	}
	int posig=0;
	for(int l=0; l<gaps.size();++l) {
		posig +=at(order[l]).size();
		mc->insertgap(posig,gaps[l]);
		posig += gaps[l];
	}
	int segposi=0;
	int l=0;
	for(int s:order) {
		int rigidstart=segposi+n_flexible_[s];
		if(l==0) rigidstart=segposi;
		int rigidend=segposi+at(s).size()-c_flexible_[s];
		if(l==order.size()-1) rigidend=segposi+at(s).size();
		mc->setrigidbetween(rigidstart,rigidend);
		if(l<gaps.size())segposi += at(s).size()+gaps[l++];
	}
}
int CoreSegments::buildclosedloops(int s1,int s2,int minlength,
		int maxlength){
	GapIndex idx=std::make_pair(s1,s2);
	if(closedloops_.find(idx) == closedloops_.end()){
		closedloops_.insert(std::make_pair(idx,LoopSolutions()));
	}
	std::vector<int> order;
	order.push_back(s1);
	order.push_back(s2);
	for(int s=0;s<size();++s) {
		if(s==s1 || s== s2) continue;
		order.push_back(s);
	}
	std::vector<int> gaps;
	gaps.resize(order.size()-1,7);
	EnergyTerms eterms;
	eterms.addstericclash(0.0);
	eterms.addtorsionene();
	EnergyTerms eterms_rank;
	eterms_rank.addstericclash(0.0);
	eterms_rank.addtorsionene();
//	eterms_rank.addtetrasefene(0.17);
	for(int len=minlength; len<=maxlength;++len) {
		gaps[0]=len;
		MainChain mc;
		buildmainchain(order,gaps,&mc);
		for( int leftext=0;leftext<c_flexible_[s1]+1; ++leftext){
			for (int rightext=0; rightext<n_flexible_[s2]+1;++rightext){
				int loopstart=sizeofsegment(s1) - leftext-1;
				int loopend=sizeofsegment(s1)+len+rightext + 1;
				ClosingLoop cl(mc, loopstart, loopend, ClosingLoop::MUTATESEQ_CONF);
				int ntotalsol=0;
				int nstop=10;
				NSPdstl::TopN<SegmentPtr> topn(nstop);
				for (int ntry = 0; ntry < 500*nstop; ++ntry) {
					int nsol = cl.solve();
					for (int i = 0; i < nsol; ++i) {
						SegmentPtr newloopi(
								new std::vector<BackBoneSite>());
						cl.getsitessolution(i, newloopi.get());
						std::map<std::string, double> energies;
						double eclash = calcloopenergy(mc, loopstart, newloopi->size(),
								*newloopi, eterms, energies, eterms.emaxallowed);

//						if(eclash >= eterms.emaxallowed || energies["TorsionEne"] >-1.0*(double) (newloopi->size())) continue;
						if(eclash >= eterms.emaxallowed) continue;
//						 checklooptorsions(mc,loopstart, newloopi->size(),*newloopi);
			//			std::cout << energies["TorsionEne"];
			//			energies.clear();
			//			double erank=calcloopenergy(mc, loopstart, newloopi->size(),
			//					*newloopi, eterms_rank, energies, eterms_rank.emaxallowed);
			//			std::cout <<"\t" <<eclash <<std::endl;
						topn.push(newloopi,eclash);
						++ntotalsol;
					} //process solution
					if(ntry>=10*nstop && ntotalsol >=nstop) break;
				} //solution tries
				std::vector<std::pair<SegmentPtr,double>> res=topN2vector(topn);
				std::cout <<"Loop energies:";
				auto & lps=closedloops_.at(idx);
				for (auto &r:res) {
					lps.push_back(LoopSolution());
					LoopSolution& ls=lps.back();
					ls.loop=r.first;
					ls.energy=r.second;
					ls.nextension=leftext;
					ls.cextension=rightext;
					std::cout <<"\t" << r.second;
				}
				std::cout << std::endl;
			} // right extension
		} //left extension
	} //loop lengths
	if(closedloops_[idx].empty()) closedloops_.erase(idx);
	return closedloops_[idx].size();
}
int CoreSegments::buildlinkedchains(const std::vector<int> & order){
	for (int i = 0; i < order.size() - 1; ++i) {
		if (closedloops_.at(std::make_pair(order[i], order[i + 1])).empty()) {
		 return 0;
		}
	}
	int nlc_old=linkedchains_.size();
	LinkedChain lc;
	lc.segments=order;
	extendchain(lc);
	return linkedchains_.size()-nlc_old;
}
void CoreSegments::extendchain(LinkedChain lc){
	if(lc.loopidx_.size()==lc.segments.size()-1){
		std::string lcname=linkedchainname(lc);
		if(builtchains_.find(lcname) != builtchains_.end()) return;
		builtchains_.insert(lcname);
		linkedchains_.push_back(lc);
		return;
	}
	StericClash *efunc = &(eneterminstance<StericClash>());
	int s1=lc.segments[lc.loopidx_.size()];
	int s2=lc.segments[lc.loopidx_.size()+1];
	GapIndex g=std::make_pair(s1,s2);
	LoopSolutions & solutions=closedloops_.at(g);
	int solidx=0;
	for(auto &sol:solutions) {
		bool accept = true;
		const std::vector<BackBoneSite> *loopcan = sol.loop.get();
		for (int sl = 0; sl < lc.loopidx_.size(); sl++) {
			int ss1 = lc.segments[sl];
			int ss2 = lc.segments[sl + 1];
			const std::vector<BackBoneSite> *loopsl =
					closedloops_.at(std::make_pair(ss1,ss2))
					[lc.loopidx_[sl]].loop.get();
			double e = looploopenergy(*loopcan, *loopsl, efunc,
					efunc->eclash());
			if (e >= efunc->eclash()) {
				accept = false;
				break;
			}
		}
		if(accept) {
			lc.loopidx_.push_back(solidx);
			double eold=lc.energy;
			lc.energy = eold+solutions[solidx].energy;
			extendchain(lc);
			lc.loopidx_.pop_back();
			lc.energy =eold;
		}
		++solidx;
	}
}
void CoreSegments::sortlinkedchains(){
	std::vector<double> ene;
	for(auto &lc:linkedchains_){
		ene.push_back(linkedchainaverageene(lc));
	}
	NSPdstl::sort12(ene,linkedchains_);
}
