/*
 * loopccd.h
 *
 *  Created on: 2016年11月3日
 *      Author: hyliu
 */

#ifndef LOOPCCD_H_
#define LOOPCCD_H_
#include "geometry/loopccd.h"
#include "geometry/crdtree.h"
#include "geometry/quatfit.h"
#include "dstl/randomengine.h"
#include <map>

namespace NSPgeometry {
template <typename ATOMKEY>
class FixEndsLoop {
public:
	typedef TopoTree<ATOMKEY> Topology;
	typedef CrdTree<ATOMKEY> Coordinate;
	FixEndsLoop(Topology *t): topo_(t){;}
	template<typename RNG>
	void setRandomRotations(Coordinate *c, RNG &rng) const {
		for(auto &r:rotatables_){
			double t=rng()*2.0*3.14159265358979323846;
			c->rotateBond(r,t);
		}
	}
	void copyStartCrdsTo(Coordinate *c) {
		c->copyCrdMap(startcrds_);
	}
	template <typename SCORER>
	bool MCclosure(Coordinate *c,int maxcycles,double rmsdtolerance,SCORER &scorer,double *rmsd=nullptr){
		std::vector<XYZ> fixpoints;
		std::vector<XYZ> movepoints;
		auto & crdmap=c->crdmap(true);
		for(auto & e:endcrds_) {
			fixpoints.push_back(e.second);
			movepoints.push_back(crdmap.at(e.first));
		}
		double dev=NSPgeometry::rmsd(fixpoints,movepoints);
		int cycles=0;
		std::map<ATOMKEY,typename Topology::Tree*> & subtrees=topo_->subtreemap();
		NSPdstl::RandomEngine<> &rneg= NSPdstl::RandomEngine<>::getinstance();
		rneg.setrealrng(0,1);
		rneg.setintrng(0,rotatables_.size()-1);
		double maxstep;
		while(dev > rmsdtolerance && cycles < maxcycles) {
			int k=rneg.intrng()();
			auto it=rotatables_.begin();
			for(int m=0;m<k; ++m) ++it;
			ATOMKEY ka=*it;
			bool acceptable=true;
			double oldscore=scorer(ka,0.0);
			maxstep=dev*0.1;
			if(maxstep > 0.1) maxstep=0.1;
			const auto & aijk=topo_->aijk(ka);
			const ATOMKEY & ja=aijk.ka;
			Line axis(crdmap[ja],crdmap[ka]-crdmap[ja]);
			double rot;
			double T=1.0;
			if(dev < 0.0) {
				rot=ccdangle(axis,fixpoints,movepoints);
				double newscore=scorer(ka,rot);
				acceptable= rneg.realrng()() <exp((oldscore-newscore)/T);
			} else {
				do {
					rot=(2.0*rneg.realrng()()-1.0)*maxstep;
//					std::cout <<"calc new score." <<std::endl;
					double newscore=scorer(ka,rot);
//					std::cout << newscore <<" " <<oldscore <<std::endl;
					acceptable= rneg.realrng()() < exp((oldscore-newscore)/T);
				}while (!acceptable);
			}
			double newrmsd=rotationtormsd(axis,rot,fixpoints,movepoints);
			double diff=exp(80*(dev-newrmsd));
//			diff=1.0;
			if( cycles == cycles/1000*1000)
				std::cout <<"cycles: " <<cycles<< "  "<< oldscore << " " <<dev <<" " <<newrmsd<< " "<<diff<<std::endl;
			if(acceptable && rneg.realrng()() < diff) {
				c->rotateBond(ka,rot);
				c->calcXYZ();
				movepoints.clear();
				for(auto & e:endcrds_) {
					movepoints.push_back(c->crdmap().at(e.first));
				}
				dev=newrmsd;
			}
			++cycles;
		}
		if(dev<=rmsdtolerance) {*rmsd=dev;return true;}
		std::cout <<"Unsuccessful closing: cycles " <<cycles << " rmsd " <<dev <<std::endl;
		return false;
	}

	bool closeEnd(Coordinate *c,int maxcycles,double rmsdtolerance,double *rmsd=nullptr){
		std::vector<XYZ> fixpoints;
		std::vector<XYZ> movepoints;
		for(auto & e:endcrds_) {
			fixpoints.push_back(e.second);
			movepoints.push_back(c->crdmap().at(e.first));
		}
/*		std::vector<XYZ> temppoints(movepoints);
		QuatFit fit;
		double rmsdfit=sqrt(fit.fitting(fixpoints,temppoints));
		std::cout <<"RMSDFIT: " <<rmsdfit <<std::endl;
		for(int m=0; m<fixpoints.size(); ++m){
			std::cout << fixpoints[m].toString() <<std::endl;
			std::cout <<temppoints[m].toString() <<std::endl;
		}
		*/
		int cycles=0;
		auto & crdmap=c->crdmap();
		double dev=NSPgeometry::rmsd(fixpoints,movepoints);
		std::map<ATOMKEY,typename Topology::Tree*> & subtrees=topo_->subtreemap();
		while(dev > rmsdtolerance && cycles < maxcycles) {
			for(auto & rt:rotatables_) {
				const auto & aijk=topo_->aijk(rt);
				const ATOMKEY & ja=aijk.ka;
				Line axis(crdmap[ja],crdmap[rt]-crdmap[ja]);
				double t=ccdangle(axis,fixpoints,movepoints);
//				std::cout <<"Rotation: " << t <<std::endl;
				c->rotateBond(rt,t);
				c->calcXYZ();
				movepoints.clear();
				for(auto & e:endcrds_) {
					movepoints.push_back(c->crdmap().at(e.first));
				}
				dev=NSPgeometry::rmsd(fixpoints,movepoints);
				if(dev <= rmsdtolerance) break;
//				std::cout <<"new rmsd: " << dev <<std::endl;
			}
			++cycles;
		}
		if(dev<=rmsdtolerance) return true;
		std::cout <<"Unsuccessful closing: cycles " <<cycles << " rmsd " <<dev <<std::endl;
		return false;
	}
	void addStartCrd(const ATOMKEY & key, const XYZ & crd){
		std::map<ATOMKEY,typename Topology::Tree*> & subtrees=topo_->subtreemap();
		auto it=subtrees.find(key);
		if(it != subtrees.end()){
			startcrds_.insert(std::make_pair(key,crd));
		}
	}
	void addEndCrd(const ATOMKEY & key, const XYZ & crd){
		std::map<ATOMKEY,typename Topology::Tree*> & subtrees=topo_->subtreemap();
		auto it=subtrees.find(key);
		if(it != subtrees.end()){
			endcrds_.insert(std::make_pair(key,crd));
		}
	}

	void addRotatable(const ATOMKEY & key) {
		std::map<ATOMKEY,typename Topology::Tree*> & subtrees=topo_->subtreemap();
			auto it=subtrees.find(key);
			if(it != subtrees.end()){
				rotatables_.insert(key);
			}
	}

protected:
	Topology * topo_;
	std::map<ATOMKEY,XYZ> startcrds_;
	std::map<ATOMKEY,XYZ> endcrds_;
	std::set<ATOMKEY> rotatables_;
};
}



#endif /* LOOPCCD_H_ */
