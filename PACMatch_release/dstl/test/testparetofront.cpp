/*
 * testparetofront.cpp
 *
 *  Created on: 2017年4月28日
 *      Author: hyliu
 */

#include "dstl/paretofront.h"
#include "dstl/randomengine.h"

#include <iostream>
#include <memory>
using namespace NSPdstl;
void normalize(std::vector<double > &dir){
	double norm=0.0;
	for(auto d:dir) norm +=d*d;
	norm=sqrt(norm);
	for(auto &d:dir) d /=norm;
}
int main() {
	int ndim=3;
	RandomEngine<> & rneg=RandomEngine<>::getinstance();
	std::vector<std::vector <double>> dirs;
	std::vector<double> bestscores;
	for (int i=0; i<10;++i) {
		std::vector<double> dir;
		for(int d=0; d<ndim; ++d) {
			dir.push_back(rneg.realrng()());
		}
		normalize(dir);
		double bestscore=0.0;
		for(int d=0; d<ndim; ++d) {
			bestscore +=-100*(d+1.0)*dir[d];
		}
		dirs.push_back(dir);
		bestscores.push_back(bestscore);
	}
	FrontSets<std::shared_ptr<int>> front(ndim);
/*	for (int i=0;i<8; ++i) {
		std::vector<double> dir;
		for(int d=0; d<3; ++d) {
			dir.push_back(rneg.realrng()());
		}
		normalize(dir);
		front.adddirection(dir,5);
	}*/
	front.generatedirections(40,5);
	std::vector<std::vector<double>> data;

	for(int i=0; i<500000; ++i) {
		data.push_back(std::vector<double>(ndim));
		for(int d=0; d<ndim; ++d) {
			data.back()[d]=rneg.realrng(-100*(d+1.0),0)();
		}
		front.save(std::shared_ptr<int>(new int(i)),data.back());
	}
	rneg.setrealrng(0,1);
	std::cout <<front.nsaved() <<std::endl;
	auto saved=front.saved();
	for (auto &s:saved) {
		std::cout <<*(s.obj_);
		for(int d=0; d<ndim;++d) {
			std::cout <<"\t"<<s.scores_[d];
		}
		for(int d=0; d<ndim;++d) {
				std::cout <<"\t"<<data[*(s.obj_)][d];
			}
		std::cout <<std::endl;
	}
	int m=0;
	for (auto & d:dirs) {
		std::cout <<m <<"\t" <<bestscores[m]<<std::endl;
		++m;
		double score;
		auto best=front.bestnalong(d,10);
		for (int i=0;i<best.size();++i) {
			for(int d=0; d<ndim;++d) {
				std::cout <<"\t"<<data[*(best[i].first)][d];
				}
			std::cout <<"\t\t\t"<<best[i].second<< std::endl;
		}
	}
}

