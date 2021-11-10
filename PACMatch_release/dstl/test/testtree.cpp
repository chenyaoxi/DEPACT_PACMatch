/*
 * testtree.cpp
 *
 *  Created on: 2016年11月4日
 *      Author: hyliu
 */

#include "dstl/tree.h"
#include <vector>
#include <string>
#include <tuple>
#include <iostream>
using namespace NSPdstl;
typedef std::tuple<int,char,int,std::string> AtomKey;
typedef TreeNode<AtomKey,std::string> AtomNode;
typedef Tree<AtomNode> AtomTree;
typedef typename AtomTree::Pointer TreePointer;

class PrintAtomKey{

public:
	void enter(AtomNode & node) const {
		int posi,rotamer;
		char res;
		std::string atomnm;
		std::tie(posi,res,rotamer,atomnm) = node.getKey();
		std::cout << posi<<res<<rotamer<<atomnm<<" "<<node.getActualNode()<<std::endl;
	}
	void leave(AtomNode & node) const {;}
	bool abort() {return false;}
};

int main(){

	TreePointer tree=AtomTree::make_tree(AtomNode(std::make_tuple(0,'G',0,"CA"),"backbone0"));
	std::vector<std::string> backbones{"backbone1","backbone2","backbone3","backbone4"};
	AtomTree *tail=tree.get();
	int posi=1;
	for(int i=0; i<backbones.size();++i){
		TreePointer child=AtomTree::make_tree(AtomNode(std::make_tuple(posi++,'G',0,"CA"),backbones[i]));
		tail=tail->addBranch(child,true);
	}

	TreePointer sc21=AtomTree::make_tree(AtomNode(std::make_tuple(2,'A',1,"CA"),"sidechain21"));
	tree->insertBranchAt(std::make_tuple(2,'G',0,"CA"),sc21);
	TreePointer sc22=AtomTree::make_tree(AtomNode(std::make_tuple(2,'E',2,"CA"),"sidechain22"));
	tree->insertBranchAtBackbone(std::make_tuple(2,'G',0,"CA"),sc22);
	TreePointer sc41=AtomTree::make_tree(AtomNode(std::make_tuple(4,'E',1,"CA"),"sidechain41"));
	tree->insertBranchAt(std::make_tuple(4,'G',0,"CA"),sc41);
	PrintAtomKey prt;
	std::cout <<"All nodes: " <<std::endl;
	tree->traverse(prt);
	std::cout <<"BackBone nodes: " <<std::endl;
	tree->traverseBackbone(prt);

	std::cout <<"Backtrace from sc22: " <<std::endl;
	bool abort{false};
	tree->findDecendant(std::make_tuple(2,'E',2,"CA"))->backtrace(prt,abort);
	std::cout <<"Backtrace from sc41 to bb2 " <<std::endl;
	AtomKey endkey(2,'G',0,"CA");
	tree->findDecendant(std::make_tuple(4,'E',1,"CA"))->backtrace(prt,abort);

	TreePointer cuttree=tree->cutBranchAt(endkey);
	std::cout <<"Remaining Tree: " <<std::endl;
	tree->traverse(prt);
	std::cout <<"CutTree: " <<std::endl;
	cuttree->traverse(prt);
	std::cout <<"CutTree tail to head: " <<std::endl;
//	cuttree->backboneTail()->backtrace(prt);

}


