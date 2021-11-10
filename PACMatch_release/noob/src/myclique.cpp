/*
 * myclique.cpp
 *
 *  Created on: 2021年7月25日
 *      Author: yxchen
 */

/*
 * Graphic match for seeding residues, replacing the open source cliquer.
 */
#include <iostream>
#include "noob/myclique.h"

using namespace MyClique;

// extend cs (size) into cs_new (size+1)
void MyClique::extendclique_map(Graph g, map<string, Clique> &csmap)
{
	vector<Vertex> vs = g.vertexes();
	if (csmap.size() == 0)
	{
		cout << "Function extendclique should use input cs with size > 0." << endl;
		exit(1);
	}
	else
		cout << "Finding cliques with size " << csmap.begin()->second.ids.size() + 1 << endl;
	int num = 0;
	map<string, Clique> csmap_new;
	for (auto it = csmap.begin(); it != csmap.end(); it++)
	{
		auto c = it->second;
		int mei = c.me_id;
		for (auto i : vs[mei].es)
		{
			if (c.ids.count(i) > 0) continue;
			bool extend = true;
			for (auto id : c.ids)
				if (vs[id].es.count(i) == 0)
				{
					extend = false;
					break;
				}
			if (extend)
			{
				string code;
				Clique c_new;
				for (auto id : c.ids)
					c_new.ids.insert(id);
				c_new.ids.insert(i);
				for (auto i : c_new.ids)
					code += "|" + to_string(i);
				if (csmap_new.count(code) > 0) continue;
				if (vs[mei].es.size() > vs[i].es.size())
					c_new.me_id = i;
				else
					c_new.me_id = mei;
				csmap_new.insert(make_pair(code, c_new));
				num++;
			} // extend
		} // each edge for ids[0]
	}
	csmap = csmap_new;
	cout << " There are " << num << " cliques." << endl;
}

// find clique in certain size
vector<Clique> MyClique::findclique_map(Graph g, int size)
{
	map<string, Clique> csmap;
	for (int i = 0; i < g.size(); i++)
	{
		Clique c;
		c.ids.insert(i);
		c.me_id = i;
		csmap.insert(make_pair(to_string(i), c));
	}
	int s_now = 1;
	while(s_now < size)
	{
		extendclique_map(g, csmap);
		if (csmap.size() == 0)
			break;
		else
			s_now++;
	}
	vector<Clique> cs;
	for (auto it = csmap.begin(); it != csmap.end(); it++)
		cs.push_back(it->second);
	return cs;
}

