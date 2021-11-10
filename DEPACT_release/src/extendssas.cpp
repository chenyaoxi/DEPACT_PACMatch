/*
 * extendssas.cpp
 *
 *  Created on: 2018年8月9日
 *      Author: hyliu
 */
#include "mmmatches.h"
extern "C" {
#include "cliquer/cliquer.h"
}
using namespace myobcode;
extern "C" {
struct UserData {
	std::vector<std::shared_ptr<SubstrAlignment>> *oldssas;
	std::vector<std::shared_ptr<SubstrAlignment>> *newssas;

};
static boolean addmergedssa(set_t s, graph_t *g, clique_options *opts) {
	std::vector<std::shared_ptr<SubstrAlignment>> &oldssas =
			*(((UserData *) opts->user_data)->oldssas);
	std::vector<std::shared_ptr<SubstrAlignment>> &ssas =
			*(((UserData *) opts->user_data)->newssas);
	int i = -1;
	std::vector<int> mmber;
	while ((i = set_return_next(s, i)) >= 0) {
		mmber.push_back(i);
	}
	std::shared_ptr<SubstrAlignment> cssa = oldssas[mmber[0]];
	for (int i = 1; i < mmber.size(); ++i) {
		cssa = std::shared_ptr < SubstrAlignment
				> (new SubstrAlignment(*cssa, *(oldssas[mmber[i]])));
	}
/*	bool redundant = false;
	for (auto &ossa : ssas) {
		if (equivalent(*cssa, *ossa)) {
			ossa = std::shared_ptr < SubstrAlignment
					> (new SubstrAlignment(*cssa, *ossa));
			redundant = true;
			break;
		}
	}
	if (!redundant) {*/
		ssas.push_back(cssa);
//	}
	return true;
}
}

void myobcode::extendssas(
		std::vector<std::shared_ptr<SubstrAlignment>> & ssas) {
	if(ssas.empty()) return;
	int size = ssas.size();
	graph_t *g = graph_new(size);
	for (int i = 0; i < size - 1; ++i) {
		for (int j = i + 1; j < size; ++j) {
			if (compatible(*(ssas[i]), *(ssas[j]))) {
				GRAPH_ADD_EDGE(g, i, j);
			}
		}
	}
	std::vector<std::shared_ptr<SubstrAlignment>> oldssas=ssas;
	UserData userdata;
	userdata.oldssas = &oldssas;
	userdata.newssas=&ssas;
	ssas.clear();
	clique_default_options->user_function = addmergedssa;
	clique_default_options->user_data = (void *) (&userdata);
	clique_default_options->output=stderr;
	clique_find_all(g, 1, size,true, NULL);
	graph_free(g);
}

