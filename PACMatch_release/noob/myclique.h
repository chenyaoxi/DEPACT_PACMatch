/*
 * myclique.h
 *
 *  Created on: 2021年7月25日
 *      Author: yxchen
 */

#ifndef NOOB_MYCLIQUE_H_
#define NOOB_MYCLIQUE_H_

#include <map>
#include <set>
#include <vector>
using namespace std;

namespace MyClique
{

struct Vertex
{
	set<int> es; // vertex that having an edge to
};
/*
struct Edge
{
	int i; // vertex_i
	int j; // vertex_j
	bool e; // true: has an edge; false: no edge
};
*/
class Graph
{
public:
	void graph_new(int vertex_num) {
		vertexes_.resize(vertex_num);
		v_num_ = vertex_num;
	}; // establish a graph
	void add_edge(int i, int j) {
		vertexes_[i].es.insert(j);
		vertexes_[j].es.insert(i);
	}; // add edge
	void graph_free() {vertexes_.clear();};
	int size() {return v_num_;};
	vector<Vertex> vertexes() {return vertexes_;};
private:
	vector<Vertex> vertexes_;
	int v_num_;
//	vector<Edge> edges_;
};

struct Clique
{
	set<int> ids; // vertex ids forming a clique.
	int me_id; // id of the min number of edges.
};

void extendclique_map(Graph g, map<string, Clique> &csmap); // idlink, Clique.
vector<Clique> findclique_map(Graph g, int size); // find clique in certain size
}



#endif /* NOOB_MYCLIQUE_H_ */
