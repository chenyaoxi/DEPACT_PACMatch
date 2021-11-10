/*
 * tree.h
 *
 *  Created on: 2016年11月4日
 *      Author: hyliu
 */

#ifndef DSTL_TREE_H_
#define DSTL_TREE_H_
#include <vector>
#include <map>
#include <memory>
#include <cassert>
namespace NSPdstl{
template <typename NODE>
class Tree {
public:
	typedef typename NODE::KeyType KeyType;
	typedef NODE NodeType;
	typedef typename NODE::ActualNodeType ActualNodeType;
	typedef std::shared_ptr<Tree> Pointer;
	static Pointer make_tree(const NODE & node){
		return Pointer(new Tree (node));
	}
	NODE & getNode(){return node_;}
	const NODE &getNode() const {return node_;}
	Tree *getParent() const {return parent_;}
	void setParent(Tree * p) {parent_=p;}
	template <typename NODEACTION>
	void setParent(Tree *p, NODEACTION & nodeaction) {
		parent_=p;
		nodeaction(node_,p);
	}
	const KeyType & getNodeKey() const {return node_.getKey();}
	KeyType & getNodeKey() {return node_.getKey();}
	Tree(const NODE & node):node_(node){;}
	Tree *addBranch(Pointer & child, bool backbone=false){
		const KeyType & childkey=child->getNodeKey();
		branches_.insert(std::make_pair(childkey,child));
		child->setParent(this);
		if(backbone) backbonekey_= &childkey;
		return child.get();
	}
	Tree *addChild(KeyType childkey, bool backbone=false){
		Pointer child=make_tree(NODE(childkey));
		return addBranch(child,backbone);
	}

	Tree *findDecendant(const KeyType & query) {
		if(getNodeKey() == query) return this;
		for(auto & bp:branches_){
			Tree *tp=bp.second->findDecendant(query);
			if(tp) return tp;
		}
		return nullptr;
	}
	Tree *findBackboneDecendant(const KeyType & query) {
		if(getNodeKey() == query) return this;
		if(!backbonekey_) return nullptr;
		if(branches_.find(*backbonekey_) != branches_.end()) {
				return branches_.at(*backbonekey_)->findBackboneDecendant(query);
		}
		return nullptr;
	}

	Tree * insertBranchAt(const KeyType &posi, Pointer child,bool backbone=false) {
		Tree *tp=findDecendant(posi);
		if(tp) {
			tp->addBranch(child,backbone);
			return child.get();
		}
		return nullptr;
	}

	Pointer cutBranchAt(const KeyType &posi) {
		Tree *tp=findDecendant(posi);
		if(tp) {
			Tree * cutposi=tp->getParent();
			Pointer cutpart=cutposi->branches_.at(posi);
			cutposi->branches_.erase(posi);
			cutpart->setParent(nullptr);
			return cutpart;
		}
		return Pointer(nullptr);
	}

	Tree * insertBranchAtBackbone(const KeyType &posi, Pointer child,bool backbone=false) {
		Tree *tp=findBackboneDecendant(posi);
		if(tp) {
			tp->addBranch(child,backbone);
			return child.get();
		}
		return nullptr;
	}

	Tree *backboneTail(){
		if(!backbonekey_) return this;
		return branches_.at(*backbonekey_)->backboneTail();
	}
	Tree *head(){
		if(!parent_) return this;
		return parent_->head();
	}
	template <typename ACTION>
	void traverse(ACTION & action){
		action.enter(node_);
		for(auto & bp:branches_) {
			bp.second->traverse(action);
		}
		action.leave(node_);
	}
	template <typename ACTION>
	void traverse(ACTION & action,bool & abort){
		action.enter(node_);
		abort=action.abort();
		if(!abort && !action.stopbranch()) {
			for(auto & bp:branches_) {
				bp.second->traverse(action,abort);
				if(abort) break;
			}
		}
		action.leave(node_);
	}

	template <typename ACTION>
	void traverseBackbone(ACTION &action) {
		action.enter(node_);
		if(backbonekey_) branches_.at(*backbonekey_)->traverseBackbone(action);
		action.leave(node_);
	}
	template <typename ACTION>
	void traverseBackbone(ACTION &action,bool & abort) {
		action.enter(node_);
		abort=action.abort();
		if(!abort)	if(backbonekey_) branches_.at(*backbonekey_)->traverseBackbone(action,abort);
		action.leave(node_);
	}


	void buildsubtreemap(std::map<KeyType,Tree<NODE> *> *smap) {
		KeyType key=getNodeKey();
		assert(smap->find(key) == smap->end());
		smap->insert(std::make_pair(key,this));
		for(auto & bp:branches_) {
			bp.second->buildsubtreemap(smap);
		}
	}

	void collectkeys(std::vector<KeyType> *keys){
		keys->push_back(getNodeKey());
		for(auto & bp:branches_) {
			bp.second->collectkeys(keys);
		}
	}

	template <typename ACTION>
	void backtrace(ACTION & action, bool &abort){
		action.enter(node_);
		abort=action.abort();
		if(!abort) if(parent_) parent_->backtrace(action,abort);
		action.leave(node_);
	}
	Tree *childBranch(const KeyType & key) {
		if(branches_.find(key) != branches_.end()) {
			return branches_[key].get();
		} else {
			return nullptr;
		}
	}
	std::map<KeyType,Pointer> & children() {return branches_;}
	const std::map<KeyType,Pointer> & children() const {return branches_;}
private:
	NODE node_;
	std::map<KeyType,Pointer> branches_;
	Tree *parent_{nullptr};
	const KeyType *backbonekey_ {nullptr};
};

template <typename KEY, typename ACTUALNODE>
class TreeNode {
public:
	typedef ACTUALNODE ActualNodeType;
	typedef KEY KeyType;
	TreeNode(const KEY & k,const ACTUALNODE & n): key_(k),actualnode_(n){;}
	TreeNode(const KEY & k):key_(k),actualnode_(k){;}
	KEY & getKey() {return key_;}
	const KEY & getKey() const {return key_;}
	ACTUALNODE & getActualNode() {return actualnode_;}
	const ACTUALNODE & getActualNode() const {return actualnode_;}
private:
	KEY key_;
	ACTUALNODE actualnode_;
};
template<typename NODE>
struct TreeMapped {
	typedef typename NODE::KeyType NodeKeyType;
	typedef std::shared_ptr<Tree<NODE>> TreePointer;
	TreeMapped(TreePointer const & t): tree(t) {buildsubtreemap();}
	TreeMapped(const NODE &rootnode){
		tree=Tree<NODE>::make_tree(rootnode);
		buildsubtreemap();
	}
	TreeMapped(const NodeKeyType &rootkey){
		tree=Tree<NODE>::make_tree(NODE(rootkey));
		buildsubtreemap();
	}
	void buildsubtreemap(){
		subtreemap.clear();
		tree->buildsubtreemap(&subtreemap);
	}

	Tree<NODE> *insertChild(const NodeKeyType & parentkey, const NodeKeyType &childkey, bool backbone=false){
		assert(hasNode(parentkey) && ! hasNode(childkey));
		Tree<NODE>* child=subtreemap[parentkey]->addChild(childkey,backbone);
		subtreemap.insert(std::make_pair(childkey,child));
		return child;
	}
	Tree<NODE> *insertBranch(const NodeKeyType & parentkey, TreePointer & child, bool backbone=false){
		assert(hasNode(parentkey));
		subtreemap[parentkey]->addBranch(child,backbone);
		child->buildsubtreemap(&subtreemap);
		return child.get();
	}

	TreePointer cutBranchAt(const NodeKeyType &posi) {
		auto it=subtreemap.at(posi);
		if(it == subtreemap.end()) return TreePointer(nullptr);
		Tree<NODE>* cut=it->second;
		std::vector<NodeKeyType> cutkeys;
		cut->collectkeys(&cutkeys);
		for(auto &k:cutkeys) subtreemap.erase(k);
		return cut->cutBranchAt(posi);
	}

	bool hasNode(const NodeKeyType & key) const { return subtreemap.find(key) != subtreemap.end();}
	Tree<NODE> * branch(const NodeKeyType &key) const {
		auto it=subtreemap.find(key);
		if(it != subtreemap.end()) return it->second;
		else return nullptr;
	}
	Tree<NODE> * branch(const NodeKeyType & key) {
		auto it=subtreemap.find(key);
		if(it != subtreemap.end()) return it->second;
		return nullptr;
	}
	std::shared_ptr<Tree<NODE>> tree{nullptr};
	std::map<NodeKeyType,Tree<NODE> *> subtreemap;
};
}



#endif /* DSTL_TREE_H_ */
