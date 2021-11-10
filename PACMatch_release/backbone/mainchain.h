/*
 * mainchain.h
 *
 *  Created on: 2017年4月22日
 *      Author: hyliu
 */

#ifndef BACKBONE_MAINCHAIN_H_
#define BACKBONE_MAINCHAIN_H_
#include "backbone/backbonesite.h"
namespace NSPproteinrep {
class Steric_enefunct {
public:
	int minsep() const {return 2;}
	double twobody(const BackBoneSite &s1, const BackBoneSite &s2) const {
		if(atomsclashed(s1,s2)) return 1.e21;
		else return 0.0;
	}
	double emaxallowed() const {return 1.e20;}
};
class MainChain: public std::vector<BackBoneSite> {
public:
	struct RigidSegments {
		std::vector<std::pair<int,int>> segments;
		void addsegment(int left,int length){
			assert(left >=0 && length>=1);
			segments.push_back(std::make_pair(left,length));
		}
		int insidesegment(int posi) const {
			int i=0;
			for(auto & s:segments) {
				if(posi > s.first && posi <s.first+s.second-1) return i;
				++i;
			}
			return -1;
		}
		void clear() {segments.clear();}
	};
	MainChain(){;}
	void generaterandom(int length);
	void insertgap(int posi, int length);
	void resetresseq();
	int maskedormissing(int posi) {
		if(missing(posi)) return 1;
		if(mask_.empty()) return 0;
		return mask_.at(posi);
	}
	void copysegment(int startposi, int length, MainChain *out, int outposi=0) const;
	void clearmask() {mask_.clear();}
	void masksegments(const std::vector<std::pair<int,int>> & segments){
		if(mask_.empty()) mask_.resize(size(),0);
		for(auto &s:segments){
			for(int i=s.first; i<s.first+s.second; ++i) mask_.at(i)=1;
		}
	}
	void umasksegments(const std::vector<std::pair<int,int>> & segments){
		if(mask_.empty()) {
			mask_.resize(size(),0);
			return;
		}
		for(auto &s:segments){
			for(int i=s.first; i<s.first+s.second; ++i) mask_.at(i)=0;
		}
	}
	void maskposition(int posi, int val=1){
		if(mask_.empty()) mask_.resize(size(),0);
		mask_.at(posi)=val;
	}
	bool missing(int posi) const {
		if(posi >= size() || posi<0) return true;
		return this->at(posi).resname=="MISSING";
	}
	std::string residuename(int posi) const {
		return this->at(posi).resname;
	}
	bool localrigid(int posi) const {
		return rigidsegments_.insidesegment(posi) >=0;
	}
	void setrigidbetween(int left, int right){
		assert(left < size() && right <=size());
		rigidsegments_.addsegment(left,right-left);
	}
	void clearrigid() { rigidsegments_.clear();}
	void write(std::ostream &os);
	void write(const std::string &filename) {
		std::ofstream ofs;
		ofs.open(filename.c_str());
		write(ofs);
		ofs.close();
	}
	bool read(std::istream &is);
	bool read(const std::string &filename) {
		std::ifstream ifs;
		ifs.open(filename.c_str());
		return read(ifs);
	}
private:
	RigidSegments rigidsegments_;
	mutable std::vector<int> mask_;
};
class ChangeRegionSelector {
public:
	ChangeRegionSelector(const MainChain *mc);
	bool selectinoneloop (int *startposi, int *endposi) const;
	bool select4seg(int *startposi, int *endposi) const;
	int flexiblebegin() const {return flexiblesites_[0];}
	int flexibleend() const {return flexiblesites_.back();}
private:
	const MainChain *mc_;
	std::vector<int> flexiblesites_;
};
void mainchainfromsegments(const std::vector<std::vector<BackBoneSite>> &segments,
		const std::vector<int> & order,const std::vector<int>& looplengths,MainChain *mc);
}

#endif /* BACKBONE_BUILDSEGMENT_H_ */
