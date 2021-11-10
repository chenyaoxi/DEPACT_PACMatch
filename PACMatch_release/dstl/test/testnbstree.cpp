/*
 * testnbstree.cpp
 *
 *  Created on: 2018年6月27日
 *      Author: hyliu
 */

#include "dstl/nbstree.h"
#include "dstl/randomengine.h"
#include <algorithm>
using namespace NSPdstl;

int main(){
	std::vector<double> startval{-5,1,20,8};
	std::vector<double> binwidth{0.1,1,5,0.2};

	NBSTree<int> tree(startval,binwidth);
	std::vector<std::vector<double>> points;
	auto & rng=RandomEngine<>::getinstance();
	for(int i=0;i<500000;++i){
		points.push_back(std::vector<double>());
		auto & p=points.back();
		for(int d=0;d<4;++d){
			double r=startval[d]+rng.realrng(0,100.0)()*binwidth[d];
			p.push_back(r);
		}
	}

	for(int i=0;i<points.size();++i){
		tree.addpoint(i,points[i]);
	}

	for(int m=0;m<10;++m){
		NBSQuery q;
		for(int d=0;d<4;++d){
			q.vals.push_back(startval[d]+rng.realrng(0,100.0)()*binwidth[d]);
			q.lcuts.push_back(2.6*binwidth[d]);
			q.rcuts.push_back(1.9*binwidth[d]);
		}
		std::vector<int> nbrs;
		tree.findneighbors(q,&nbrs);
		std::sort(nbrs.begin(),nbrs.end());
		for(auto i:nbrs){
			std::cout << " "<<i;
		}
		std::cout<<std::endl;
		for(int i=0;i<points.size();++i){
			bool keep=true;
			for(int d=0;d<4;++d){
				double delta=points[i][d]-q.vals[d];
				if(-(q.lcuts[d]) <=delta && delta<= q.rcuts[d])
					continue;
				else {
					keep=false;
					break;
				}
			}
			if(keep) std::cout<<" "<<i;
		}
		std::cout<<std::endl;
	}

}


