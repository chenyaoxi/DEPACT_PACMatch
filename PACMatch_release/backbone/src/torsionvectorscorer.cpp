/*
 * torsionvectorscorer.cpp
 *
 *  Created on: 2016年12月28日
 *      Author: hyliu
 */

#include "backbone/torsionvectorscorer.h"

using namespace NSPproteinrep;
using namespace domaintree;
const double TorsionVectorScorer::motifunfoundscore{0};
const double TorsionVectorScorer::d2max{625.0};
const double TorsionVectorScorer::d20{25.0};
const double TorsionVectorScorer::resolution{2.0};
void TorsionVectorScorer::init(std::vector<BackBoneSite> *sites,int length) {
	templatevectors_= std::shared_ptr<std::vector<std::vector<double>>> (new std::vector<std::vector<double>>());
	length_=length;
	for(auto bsiter=sites->begin()+1; bsiter !=sites->end(); ++bsiter){
		if(!fragstartsite(bsiter,sites->end(),length)) continue;
		bool containcoil=false;
		for(int i=0;i<length;++i) {
			auto it=bsiter+i;
			if(it->sscodechar() == 'C') {
				containcoil=true;
				break;
			}
		}
		if(!containcoil) continue;
		templatevectors_->push_back(std::vector<double>());
		std::vector<double> &tv=templatevectors_->back();
		for(int i=0;i<length; ++i) {
			auto it=bsiter+i;
			tv.push_back(it->phi());
			tv.push_back(it->psi());
		}
		std::string motifname=getmotifname(bsiter,length);
		if(templatetrees_.find(motifname) == templatetrees_.end()) {
			templatetrees_.insert(std::make_pair(motifname,TorsionVectorTree()));
			TorsionVectorTree & tree=templatetrees_.at(motifname);
			tree.init(templatevectors_,resolution);
		}
		TorsionVectorTree & tree=templatetrees_.at(motifname);
		auto it=bsiter;
		long location=templatevectors_->size()-1;
		tree.insertpoint(location);
	}
}
double TorsionVectorScorer::score(const std::vector<double> & torsions, const std::vector<double> &tmpl) {
	double s=1.0;
	for (int i=0;i<2*length_;++i) {
		double si;
		double diff2=torsions[i]-tmpl[i];
		while(diff2 <-180.0) diff2+=360.0;
		while(diff2 >180.0) diff2 -=360.0;
		diff2=diff2*diff2;
		if(diff2 <= d20) si=1.0;
		else if(diff2 >=d2max) si=0.0;
		else {
//			double x=(diff2-d20)/(d2max-d20);
//			si=1-x*x;
			double x=(diff2-d2max)/(d20-d2max);
			si=x*x;
		}
		s*=si;
	}
	return s;
}
double TorsionVectorScorer::score(const std::string &motifname, const std::vector<double> & torsions){
	if(templatetrees_.find(motifname) == templatetrees_.end()) return motifunfoundscore;
	typename TorsionVectorTree::Tree &tree=templatetrees_.at(motifname).gettree();
	double bound2=neighborcut2();
	domaintree::D2Leaf<long,std::vector<std::vector<double>>,AngleCrd>
		d2leaf(templatevectors_.get(),1000000,bound2);
/*	for (unsigned int i=0; i<torsions.size(); ++i){
		std::cout <<torsions[i] <<"\t";
	}
	std::cout <<std::endl;*/
	tree.findneighbor(torsions,d2leaf,bound2);
	std::vector<std::pair<long,double>> &neighbors=d2leaf.nnearest().neighbors();
	double s=motifunfoundscore;
	double alpha=1.e-7;
	for(auto &n:neighbors){
		if(n.second < 0.00001) continue;  //ignore self
		s+=score(torsions,(*templatevectors_)[n.first]);
	}
	s = alpha+s/templatetrees_.at(motifname).size();
	return -log(s)-9.3;
}

double TorsionVectorScorer::score(std::vector<BackBoneSite>::iterator iter) {
	std::vector<double> tv;
	std::string motifname=getmotifname(iter,length_);
	for(unsigned int i=0;i<length_; ++i){
		tv.push_back(iter->phi());
		tv.push_back(iter->psi());
		++iter;
	}
	return score(motifname,tv);
}

double TorsionVectorScorer::scorerange(std::vector<BackBoneSite>::iterator begin, std::vector<BackBoneSite>::iterator end){
	double stot=0.0;
	auto be=begin;
/*	std::vector<double> count(end-begin,0.0);
	while(end-be >=length_) {
		for(unsigned int i=0; i<length_;++i)
		   count[be-begin+i]+=1.0;
		++be;
	}
	be=begin;
	std::vector<double> w(end-begin-length_+1,0.0);
	while(end-be >=length_){
		for (unsigned int i=0; i<length_; ++i)
			w[be-begin] +=count[be-begin+i];
		w[be-begin] /=(double) length_;
		++be;
	}
	be=begin;
	*/
	while (end-be >=length_) {
		stot += score(be)/(double)(length_);
		++be;
	}
	return stot;
}


