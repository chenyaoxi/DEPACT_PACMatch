/*
 * graph.cpp
 *
 *  Created on: 2017年8月20日
 *      Author: hyliu
 */
#include "dstl/graph.h"
using namespace NSPdstl;
void MCPAlgorithm::run(const Graph &G, std::vector<int> *path){
	init(G,path);
	for(int i=gsize_-1;i>=0;--i){
		std::vector<int> inbrs;
		for(int j=i+1; j<gsize_;++j){
			if(G.isneighbor(i,j)) inbrs.push_back(j);
			dfs(inbrs,1);
		}
		bestat_[i]=bestn_;
	}
}
bool MCPAlgorithm::dfs(const std::vector<int> & nbrs, int cnt){
	if(nbrs.empty()) {
		if(bestn_<cnt){
			bestn_=cnt;
			path_->clear();
			for(int i=0;i<cnt;++i){
				path_->push_back(tmp_[i]);
			}
			return true;
		} else
			return false;
	}
	for(int i=0;i<nbrs.size();++i){
		if(cnt+nbrs.size()-i<=bestn_)return false;
		int vi=nbrs[i];
		if(cnt+bestat_[vi]<=bestn_) return false;
		tmp_[cnt]=vi;
		std::vector<int> inbrs;
		for(int j=i+1;j<nbrs.size();++j){
			if(g->isneighbor(vi,nbrs[j]))
				inbrs.push_back(nbrs[j]);
		}
		if (dfs(inbrs,cnt+1)) return true;
	}
	return false;
}


