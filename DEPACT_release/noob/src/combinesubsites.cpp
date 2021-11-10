/*
 * combinesubsites.cpp
 *
 *  Created on: 2018年8月13日
 *      Author: hyliu
 */

#include "subsite.h"
extern "C" {
#include "cliquer/cliquer.h"
}
using namespace subsitedesign;
using namespace NSPproteinrep;
extern "C" {
struct psubsites {
	const std::vector<std::shared_ptr<Subsite>> *oldssites;
	std::vector<std::shared_ptr<Subsite>> *newssites;
	std::vector<std::vector<int>> *mmbers;

};
static boolean addmergedssites(set_t s, graph_t *g, clique_options *opts) {
	const std::vector<std::shared_ptr<Subsite>> &oldssites =
			*(((psubsites *) opts->user_data)->oldssites);
	std::vector<std::shared_ptr<Subsite>> &ssites =
			*(((psubsites *) opts->user_data)->newssites);
	int i = -1;
	std::vector<int> mmber;
	while ((i = set_return_next(s, i)) >= 0) {
		mmber.push_back(i);
	}
	std::shared_ptr<Subsite> css = oldssites.at(mmber[0]);
	for (int i = 1; i < mmber.size(); ++i) {
		css = std::shared_ptr < Subsite
				> (new Subsite(*css, *(oldssites.at(mmber[i]))));
	}
	ssites.push_back(css);
	return true;
}
static boolean addmergedssites_mm(set_t s, graph_t *g, clique_options *opts) {
	const std::vector<std::shared_ptr<Subsite>> &oldssites =
			*(((psubsites *) opts->user_data)->oldssites);
	std::vector<std::shared_ptr<Subsite>> &ssites =
			*(((psubsites *) opts->user_data)->newssites);
	std::vector<std::vector<int>> &mmbers =
			*(((psubsites *) opts->user_data)->mmbers);
	int i = -1;
	std::vector<int> mmber;
	while ((i = set_return_next(s, i)) >= 0) {
		mmber.push_back(i);
	}
	std::shared_ptr<Subsite> css = oldssites.at(mmber[0]);
	for (int i = 1; i < mmber.size(); ++i) {
		css = std::shared_ptr < Subsite
				> (new Subsite(*css, *(oldssites.at(mmber[i]))));
	}
	ssites.push_back(css);
	mmbers.push_back(mmber);
	return true;
}
}


std::vector<std::shared_ptr<Subsite>> subsitedesign::combinebycliques(
		const std::vector<std::shared_ptr<Subsite>> & ssites){
	if(ssites.empty()) return std::vector<std::shared_ptr<Subsite>>();
	int size = ssites.size();
	graph_t *g = graph_new(size);
	for (int i = 0; i < size - 1; ++i) {
		for (int j = i + 1; j < size; ++j) {
			if (compatible(*(ssites[i]), *(ssites[j]))) {
				GRAPH_ADD_EDGE(g, i, j);
			}
		}
	}
	std::vector<std::shared_ptr<Subsite>> newssites;
	psubsites userdata;
	userdata.oldssites = &ssites;
	userdata.newssites=&newssites;
	clique_default_options->user_function = addmergedssites;
	clique_default_options->user_data = (void *) (&userdata);
	clique_default_options->output=stderr;
	clique_find_all(g, 1, size,true, NULL);
	graph_free(g);
	return newssites;
}

std::vector<std::shared_ptr<Subsite>> subsitedesign::combinebycliques_mm(
		const std::vector<std::shared_ptr<Subsite>> & ssites,
		std::vector<std::vector<int>> & mmbers){
	if(ssites.empty()) return std::vector<std::shared_ptr<Subsite>>();
	int size = ssites.size();
	graph_t *g = graph_new(size);
	for (int i = 0; i < size - 1; ++i) {
		for (int j = i + 1; j < size; ++j) {
			if (compatible(*(ssites[i]), *(ssites[j]))) {
				GRAPH_ADD_EDGE(g, i, j);
			}
		}
	}
	std::vector<std::shared_ptr<Subsite>> newssites;
	psubsites userdata;
	userdata.oldssites = &ssites;
	userdata.newssites = &newssites;
	userdata.mmbers = &mmbers;
	clique_default_options->user_function = addmergedssites_mm;
	clique_default_options->user_data = (void *) (&userdata);
	clique_default_options->output=stderr;
	clique_find_all(g, 1, size,true, NULL);
	graph_free(g);
	return newssites;
}

bool subsitedesign::compatible(const Subsite &s1, const Subsite &s2){
// clash distances squared between two atoms of different nseparatingbonds
	static const std::vector<std::vector<double>> MINDIST2{
		{12.25,12.25,12.25,4.0},
		{12.25,12.25,12.25,4.0},
		{12.25,12.25,12.25,4.0},
		{4.0,4.0,4.0,4.0}
	};
	for(auto &kvec1:s1.keyresidues()){
		for (auto &kr1:kvec1){
			const AAConformer &conf1=kr1.conformer;
			const std::vector<int> & nsep1=kr1.nseparatingbonds;
			for(int i=0;i<conf1.atomlist.size();++i){
				NSPgeometry::XYZ xi=conf1.globalcrd.at(conf1.atomlist.at(i));
				int nsepi=nsep1.at(i);
				for(auto &kvec2:s2.keyresidues()){
					for(auto &kr2:kvec2){
						const AAConformer &conf2=kr2.conformer;
						const std::vector<int> &nsep2=kr2.nseparatingbonds;
						for(int j=0;j<conf2.atomlist.size();++j){
							int nsepj=nsep2.at(j);
							double r2=(xi-conf2.globalcrd.at(conf2.atomlist.at(j))).squarednorm();
							if(r2<MINDIST2[nsepi][nsepj]) {
								return false;
							}
						} //conf2 atoms
					}//kr2
				} //kvec2
			} //conf1 atoms
		} //kr1
	}//kvec1
	return true;
}

