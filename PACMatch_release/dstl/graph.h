/*
 * graph.h
 *
 *  Created on: 2017年8月20日
 *      Author: hyliu
 */

#ifndef DSTL_GRAPH_H_
#define DSTL_GRAPH_H_
#include <vector>
namespace NSPdstl{
class Graph{
public:
	typedef int Vertex;
	Graph(int s=0):size(s){;}
	virtual bool isneighbor(Vertex v1,Vertex v2) const=0;
	virtual ~Graph(){;}
	int size;
private:

};

class MCPAlgorithm{
public:
	void run(const Graph &G,std::vector<int> *path);
private:
	bool dfs(const std::vector<int> & nbrs,int cnt);
	int bestn_;
	std::vector<int> bestat_;
	std::vector<int> *path_;
	std::vector<int> tmp_;
	const Graph *g;
	int gsize_;
	void init(const Graph &G  ,std::vector<int> *path){
		gsize_=G.size;
		g=&G;
		bestn_=0;
		path_=path;
		path_->clear();
		bestat_.clear();
		bestat_.resize(G.size,0);
		tmp_.resize(G.size,0);
	}
};

inline std::vector<int> maxclique(Graph &G){
	std::vector<int> path;
	MCPAlgorithm mcp;
	mcp.run(G,&path);
	return path;
}

}
#endif /* DSTL_GRAPH_H_ */
