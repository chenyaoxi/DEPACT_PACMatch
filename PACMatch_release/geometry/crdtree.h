/*
 * coordinatetree.h
 *
 *  Created on: 2016年11月12日
 *      Author: hyliu
 */

#ifndef GEOMETRY_CRDTREE_H_
#define GEOMETRY_CRDTREE_H_
#include "geometry/calculators.h"
#include "dstl/tree.h"
#include <map>
#include <cassert>
#include <set>
namespace NSPgeometry {


template <typename ATOMKEY>
class TopoTree{
public:
	struct AIJK;
	typedef NSPdstl::TreeNode<ATOMKEY,AIJK> TreeNode;
	typedef NSPdstl::Tree<TreeNode> Tree;
	struct AIJK{
		AIJK(const ATOMKEY &i,const ATOMKEY &j, const ATOMKEY &k): ia(i),ja(j),ka(k){;}
		AIJK(const ATOMKEY &l):ia(l),ja(l),ka(l){;}
		ATOMKEY ia;
		ATOMKEY ja;
		ATOMKEY ka;
		bool appropriate() {return ia != ja;}
		void setParent (const TreeNode & pnode) {
			ka=pnode.getKey();
			ja=pnode.getActualNode().ka;
			ia=pnode.getActualNode().ja;
		}
	};
	TopoTree(const ATOMKEY &key0):treemapped_(key0){;}

	Tree *attachAtom(const ATOMKEY & parentkey, const ATOMKEY &childkey, bool backbone=false){
		Tree *child=treemapped_.insertChild(parentkey,childkey,backbone);
		child->getNode().getActualNode().setParent(treemapped_.branch(parentkey)->getNode());
		return child;
	}

	Tree *attachBranch(const ATOMKEY & parentkey, typename Tree::Pointer &child, bool backbone=false){
		return treemapped_.insertBranch(parentkey,child,backbone);
		std::vector<ATOMKEY> decendants;
		child->collectkeys(&decendants);
		for(auto &d:decendants) {
			TreeNode &p=treemapped_.branch(d)->getParent()->getNode();
			aijk(d).setParent(p);
		}
		return child;
	}

	Tree *tree() {return treemapped_.tree.get();}
	const Tree *tree() const {return treemapped_.tree.get();}
	std::map<ATOMKEY,Tree*> & subtreemap() {return treemapped_.subtreemap;}
	const std::map<ATOMKEY,Tree*> & subtreemap() const {return treemapped_.subtreemap;}
	NSPdstl::TreeMapped<TreeNode> & treemapped() {return treemapped_;}
	const NSPdstl::TreeMapped<TreeNode> & treemapped() const {return treemapped_;}
	AIJK & aijk(ATOMKEY key) {
		return treemapped_.branch(key)->getNode().getActualNode();
	}
	const AIJK & aijk(ATOMKEY key) const {
		return treemapped_.branch(key)->getNode().getActualNode();
	}
	Tree *branch(ATOMKEY key) {
		return treemapped_.branch(key);
	}
	Tree *branch(ATOMKEY key) const {
		return treemapped_.branch(key);
	}
protected:
	NSPdstl::TreeMapped<TreeNode> treemapped_;
};


template <typename ATOMKEY> class CalcNodeInternal;
template <typename ATOMKEY> class CalcNodeXYZ;
struct IntCrd{
		IntCrd(double d=-1, double a=1000.0, double t=1000.0): distance(d),angle(a),torsion(t){;}
		double distance;
		double angle;
		double torsion;
		static bool validAngle(double a) { return a < 100;}
		static bool validTorsion(double t) {return t<100;}
		static bool validDistance(double d) {return d>0;}
		bool valid() const {return validDistance(distance) &&
				validAngle(angle) && validTorsion(torsion);
		}
		template <typename ATOMKEY>
		IntCrd(const TopoTree<ATOMKEY> * topo_,const std::map<ATOMKEY,XYZ> & crdref,const ATOMKEY &key):
		distance(-1.0),angle(10000),torsion(10000){
			const auto & itl=crdref.find(key);
			if(itl != crdref.end()) {
				assert(itl->second.valid());
				const typename TopoTree<ATOMKEY>::AIJK & aijk=topo_->aijk(key);
				if(aijk.ka != key){
					const auto & itk=crdref.find(aijk.ka);
						if(itk != crdref.end()) {
							assert(itk->second.valid());
							distance=NSPgeometry::distance(itl->second, itk->second);
							if(aijk.ja != aijk.ka) {
								const auto & itj=crdref.find(aijk.ja);
								if(itj != crdref.end()) {
									assert(itj->second.valid());
									angle=NSPgeometry::angle(itj->second, itk->second, itl->second);
									if(aijk.ia != aijk.ja) {
										const auto & iti=crdref.find(aijk.ia);
										if(iti != crdref.end()) {
											assert(iti->second.valid());
											torsion=NSPgeometry::torsion(iti->second,
											itj->second,itk->second,itl->second);
										} //iti
									} //ia
								} //itj
							} //ja
						}//itk
				} //ka
			} //itl
		}
};

template <typename ATOMKEY>
class CrdTree {
public:
	CrdTree(TopoTree<ATOMKEY> *t):topotree_(t){resetMaps();}
	std::map<ATOMKEY,IntCrd> & intcrdmap(bool tryupdate=false){
		if(tryupdate)
			if(crdupdatepoints_.empty()) calcInternal();
		return intcrdmap_;}
	const std::map<ATOMKEY,IntCrd> & intcrdmap() const {return intcrdmap_;}
	std::map<ATOMKEY,XYZ> & crdmap(bool tryupdate=false){
		if(tryupdate)
			if(!crdupdatepoints_.empty()) calcXYZ();
		return crdmap_;}
	const std::map<ATOMKEY,XYZ> & crdmap() const {return crdmap_;}
	TopoTree<ATOMKEY>* & topotree() {return topotree_;}
	TopoTree<ATOMKEY>* const & topotree() const {return topotree_;}

	void resetCrdMap(){
		crdmap_.clear();
		crdupdatepoints_.clear();
		crdupdatepoints_.insert(topotree_->tree()->getNodeKey());
		std::map<ATOMKEY,typename TopoTree<ATOMKEY>::Tree *> & subtreemap= topotree_->subtreemap();
		for(auto & st:subtreemap) {
			ATOMKEY key=st.first;
			crdmap_.insert(std::make_pair(key,XYZ(1.e20,1.e20,1.e20)));
		}
	}
	bool crdValid(const ATOMKEY & key) {
		auto it=crdmap_.find(key);
		if( it !=crdmap_.end())	return (it->second.x_< 1.e19 && it->second.y_ <1.e19 && it->second.z_<1.e19);
		return false;
	}
	void resetIntCrdMap(){
		intcrdmap_.clear();
		std::map<ATOMKEY,typename TopoTree<ATOMKEY>::Tree *> & subtreemap= topotree_->subtreemap();
		for(auto & st:subtreemap) {
			ATOMKEY key=st.first;
			intcrdmap_.insert(std::make_pair(key, IntCrd()));
		}
	}
	void resetMaps() {
		resetCrdMap();
		resetIntCrdMap();
	}
	void calcXYZ(bool forceupdate=false){
		assert (crdmap_.size() >= topotree_->subtreemap().size());
		atomsmoved_.clear();
		CalcNodeXYZ<ATOMKEY> calcnodexyz(this,forceupdate);
		topotree_->tree()->traverse(calcnodexyz);
		crdupdatepoints_.clear();
	}
	void calcInternal(){
		assert (intcrdmap_.size() >= topotree_->subtreemap().size());
		CalcNodeInternal<ATOMKEY> calcnodeinternal(this);
		topotree_->tree()->traverse(calcnodeinternal);
	}

	void setDistance(const ATOMKEY &la, double d){
		if(d != intcrdmap_[la].distance) {
			crdupdatepoints_.insert(la);
			intcrdmap_[la].distance=d;
		}
	}

	void setAngle(const ATOMKEY &la, double a){
		if(a != intcrdmap_[la].angle ) {
			crdupdatepoints_.insert(la);
			intcrdmap_[la].angle=a;
		}
	}
	void setTorsion(const ATOMKEY &la, double t){
		if( t != intcrdmap_[la].torsion) {
			crdupdatepoints_.insert(la);
			double tpi=2.*3.14159265358979323846;
			if( t>tpi ) t -=tpi;
			else if(t<-tpi) t+=tpi;
			intcrdmap_[la].torsion=t;
		}
	}
	void rotateBond(const ATOMKEY &ka,double dt){
		auto kachildren=topotree_->branch(ka)->children();
		for( auto & l:kachildren) {
			setTorsion(l.first, intcrdmap_[l.first].torsion+dt);
		}
	}
	void setInternal(const ATOMKEY &la, double d,double a,double t){
		IntCrd ic(d,a,t);
		if(intcrdmap_[la] != ic){
			crdupdatepoints_.insert(la);
			intcrdmap_[la]=IntCrd(d,a,t);
		}
	}
	void updateInternal(const std::map<ATOMKEY,XYZ> & crdref) {
		for(auto &ic:intcrdmap_){
			IntCrd nic(topotree_,crdref,ic.first);
			if(IntCrd::validDistance(nic.distance)) {
				crdupdatepoints_.insert(ic.first);
				ic.second.distance=nic.distance;
				if(IntCrd::validAngle(nic.angle)) ic.second.angle=nic.angle;
				if(IntCrd::validTorsion(nic.torsion)) ic.second.torsion=nic.torsion;
			}
		}
	}
	void addAtom(const ATOMKEY &la, const IntCrd & ic) {
		auto it=intcrdmap_.find(la);
		if( it != intcrdmap_.end()) it->second=ic;
		else {
			intcrdmap_.insert(std::make_pair(la,ic));
			crdmap_.insert(std::make_pair(la,XYZ()));
		}
		crdupdatepoints_.insert(la);
	}

	bool independentXYZ(const ATOMKEY & la){
		typename TopoTree<ATOMKEY>::AIJK & aijk=topotree_->aijk(la);
		return !aijk.appropriate();
	}

	void copyCrdMap(const std::map<ATOMKEY,XYZ> & crd0){
		updateInternal(crd0);
		for(auto & c0:crd0) {
			auto it=crdmap_.find(c0.first);
			if(it != crdmap_.end()) {
				if(crdupdatepoints_.find(c0.first) != crdupdatepoints_.end())
					crdupdatepoints_.erase(c0.first);
				it->second = c0.second;
			}
		}
	}
	bool startupdate(ATOMKEY key) {return crdupdatepoints_.find(key) != crdupdatepoints_.end();}
	std::set<ATOMKEY> & atomsmoved() {return atomsmoved_;}
	const std::set<ATOMKEY> &atomsmoved() const {return atomsmoved_;}
protected:
	TopoTree<ATOMKEY>* topotree_;
	std::map<ATOMKEY, XYZ> crdmap_;
	std::map<ATOMKEY,IntCrd> intcrdmap_;
	std::set<ATOMKEY> crdupdatepoints_;
	std::set<ATOMKEY> atomsmoved_;
};

template<typename ATOMKEY>
class CalcNodeXYZ {
public:
	typedef ATOMKEY KeyType;
	typedef typename TopoTree<KeyType> ::TreeNode Node;
	typedef typename  TopoTree<KeyType>::Tree Tree;
	CalcNodeXYZ( CrdTree<KeyType> *ct,bool forceupdate=false):crdtree_(ct),update_(forceupdate){;}
	void enter(Node & node) {
		ATOMKEY &key=node.getKey();
		if(!update_){
			if(crdtree_->startupdate(key)) {
				updatestartkey_=key;
				update_=true;
			}
		}
		if(update_){
			typename TopoTree<KeyType>::AIJK & aijk=node.getActualNode();
			if(aijk.appropriate()){
				IntCrd & ic=crdtree_->intcrdmap()[key];
				if(!ic.valid()) {
					std::cout <<ic.distance <<" " <<ic.angle <<" " <<ic.torsion <<std::endl;
				}
				assert(ic.valid());
				std::map<ATOMKEY, XYZ> & crdmap=crdtree_->crdmap();
				assert(crdmap[aijk.ka].valid() && crdmap[aijk.ja].valid() && crdmap[aijk.ia].valid());
				crdmap[key]=InternaltoXYZ(crdmap[aijk.ka],crdmap[aijk.ja],crdmap[aijk.ia],ic.distance,
						ic.angle,ic.torsion);
				crdtree_->atomsmoved().insert(key);
			}
		}
	}
	void leave(Node & node){
		ATOMKEY & key=node.getKey();
		if(key == updatestartkey_) update_=false;
		;}
private:
	CrdTree<ATOMKEY> *crdtree_;
	bool update_;
	ATOMKEY updatestartkey_;
};

template<typename ATOMKEY>
class CalcNodeInternal {
public:
	typedef ATOMKEY KeyType;
	typedef typename TopoTree<KeyType>::TreeNode  Node;
	typedef typename TopoTree<KeyType>::Tree Tree;
	CalcNodeInternal(CrdTree<KeyType> *ct):crdtree_(ct){;}
	void enter(Node & node) {
		ATOMKEY &key=node.getKey();
		typename TopoTree<KeyType>::AIJK & aijk=node.getActualNode();
		if(aijk.ka==key) return;
		IntCrd & ic=crdtree_->intcrdmap()[key];
		std::map<ATOMKEY, XYZ> & crdmap=crdtree_->crdmap();
		assert(crdmap[key].valid() && crdmap[aijk.ka].valid());
		ic.distance=distance(crdmap[aijk.ka],crdmap[key]);
		if(aijk.ka ==  aijk.ja) return;
		assert(crdmap[aijk.ja].valid());
		ic.angle=angle(crdmap[aijk.ja],crdmap[aijk.ka],crdmap[key]);
		if(aijk.ja == aijk.ia) return;
		assert(crdmap[aijk.ia].valid());
		ic.torsion=torsion(crdmap[aijk.ia],crdmap[aijk.ja],crdmap[aijk.ka],crdmap[key]);
	}
	void leave(Node & node){;}
private:
	CrdTree<KeyType> *crdtree_;
};

}


#endif /* GEOMETRY_CRDTREE_H_ */
