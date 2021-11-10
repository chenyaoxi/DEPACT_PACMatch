/*
 * paretofront.h
 *
 *  Created on: 2017年4月27日
 *      Author: hyliu
 */

#ifndef DSTL_PARETOFRONT_H_
#define DSTL_PARETOFRONT_H_
#include "topn.h"
#include "geometry/spherepoints.h"
#include <set>
#include <cstdlib>
#include <cmath>
#include <cassert>
namespace NSPdstl {
template<typename T>
class FrontSets {
public:
	struct MultiScoredType {
		T obj_;
		std::vector<double> scores_;
		MultiScoredType(){;}
		MultiScoredType(const T &obj,const std::vector<double> scores):
			obj_(obj),scores_(scores){;}
		double combscore(const std::vector<double> &dir) const {
			double res=0.0;
			for(int d=0; d<scores_.size();++d) {
				res +=scores_[d]*dir[d];
			}
			return res;
		}
		bool operator<(const MultiScoredType & st2) const {
			return obj_ <st2.obj_;
		}
	};
	struct FrontSet{
		std::vector<double> direction;
		TopN<T> topn;
		FrontSet() {;}
		FrontSet(const std::vector<double> & dir, int n):
			direction(dir),topn(n){;}
		double combscore(const std::vector<double> &scores) const {
			double res=0.0;
			for(int d=0; d<scores.size();++d) {
				res +=scores[d]*direction[d];
			}
			return res;
		}
	};
	bool save(const T & obj, const std::vector<double> &score){
		assert(score.size()== ndim_score_);
		assert(sets_.size()>0);
		bool keep=false;
		for(auto &s:sets_) {
			T toremove=obj;
			double cs=s.combscore(score);
			bool keeps=s.topn.push(obj,cs,&toremove);
			if(toremove !=obj)
			{
				MultiScoredType mt(toremove,score);
				auto it=saved_.find(mt);
				saved_.erase(it);
			}
			if(keeps) {
				saved_.insert(MultiScoredType(obj,score));
			}
			keep=keep||keeps;
		}
		return keep;
	}
	bool willsave(const std::vector<double> &score){
		assert(score.size()== ndim_score_);
		assert(sets_.size()>0);
		bool keep=false;
		for(auto &s:sets_) {
			double cs=s.combscore(score);
			bool keeps=s.topn.keep(cs);
			keep=keep||keeps;
		}
		return keep;
	}
	bool willsavetodir(const std::vector<double> &score,FrontSet *&fs,double *cs){
		fs=&(findset(score,cs));
		return fs->topn.keep(*cs);
	}
	void savetodir(const T & obj,const std::vector<double> &score,FrontSet &fs,double cs){
		T toremove=obj;
		fs.topn.push(obj,cs,&toremove);
		saved_.insert(MultiScoredType(obj,score));
		if(toremove !=obj) {
			MultiScoredType mt(toremove,score);
			auto it=saved_.find(mt);
			saved_.erase(it);
		}
	}
	FrontSets(int score_dim=2):ndim_score_(score_dim){;}
	std::set<MultiScoredType> saved() const {
		std::set<MultiScoredType> res;
		for(auto &s:saved_) res.insert(s);
		return res;
	};
	T bestalong(const std::vector<double> & dir, double *comscore) const {
		double bestscore=100000000;
		const T *bestobj;
		auto savedobj=saved();
		for(auto & s:savedobj) {
			double tmp=s.combscore(dir);
			if(tmp <bestscore) {
				bestscore=tmp;
				bestobj=&(s.obj_);
			}
		}
		*comscore=bestscore;
		return *bestobj;
	}
	std::vector<std::pair<T,double>>  bestnalong(const std::vector<double> & dir, int n) {
		TopN<T> tops(n);
		auto savedobj=saved();
		for(auto & s:savedobj) {
			double tmp=s.combscore(dir);
			tops.push(s.obj_,tmp);
		}
		return topN2vector(tops);
	}
	void adddirection(const std::vector<double> &dir, int ntop){
		sets_.push_back(FrontSet(dir,ntop));
		ndirections_=sets_.size();
	}
	int generatedirections(int ndir,int ntop) {
		if(ndim_score_==2) {
			return generatedirections2D(ndir,ntop);
		} else if (ndim_score_==3) {
			return generatedirections3D(ndir,ntop);
		} else {
			abort();
		}
	}
	int generatedirections2D(int ndir,int ntop,double thetamin=0,double thetamax=90.0) {
		double deg=3.14159265/180.0;
		double step=(thetamax-thetamin)/(double)(ndir-1);
		double theta=thetamin;
		for (int i=0; i<ndir; ++i) {
			std::vector<double> dir;
			dir.push_back(cos(theta*deg));
			dir.push_back(sin(theta*deg));
			adddirection(dir,ntop);
			theta +=step;
		}
		return sets_.size();
	}
	int generatedirections3D(int ndir,int ntop) {
		adddirection(std::vector<double>({1,0,0}),ntop);
		adddirection(std::vector<double>({0,1,0}),ntop);
		adddirection(std::vector<double>({0,0,1}),ntop);
		int npoints=8*(ndir-3);
		std::vector<std::vector<double>> dirs;
		NSPgeometry::genspherepoints(npoints, dirs);
		for (int i=0; i<npoints; ++i) {
			std::vector<double> &dir=dirs[i];
			bool skip=(dir[0] <0 || dir[1] <0|| dir[2] <0);
			if(!skip)adddirection(dir,ntop);
		}
		return sets_.size();
	}
//	const TopN<T> & getset(const std::vector<double> &dir) const {
//		return findset(dir,nullptr).topn;
//	}
	int ndirections() const {return ndirections_;}
	int nsaved() const{return saved().size();}
	const MultiScoredType &getsaved(int i) const {
		assert(i<saved_.size());
		auto iter=saved_.begin();
		for (int m=0;m<i;++m ) ++iter;
		return *iter;
	}
private:
	std::vector<FrontSet> sets_;
	std::multiset<MultiScoredType> saved_;
	int ndirections_{0};
	int ndim_score_;
	FrontSet & findset(const std::vector<double> & score,double *combscore){
		int match=-1;
		double bestmatch=-10000000;
		for (int i=0;i<ndirections_;++i) {
			double tmpmatch=sets_[i].combscore(score);
			if(tmpmatch <0) tmpmatch=-tmpmatch;
			if(tmpmatch>bestmatch) {
				bestmatch=tmpmatch;
				match=i;
			}
		}
		if(combscore){
			*combscore=sets_[match].combscore(score);
		}
		return sets_[match];
	}
};
}



#endif /* DSTL_PARETOFRONT_H_ */
