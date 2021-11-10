/*
 * exam.cpp
 *
 *  Created on: 2018年2月26日
 *      Author: hyliu
 */

#include "dataio/inputlines.h"
#include <iostream>
#include <fstream>
#include <map>
using namespace NSPdataio;

int main(int argc, char** argv){
	InputLines il;
	std::map<std::string,int> nscores;
	std::map<std::string,std::vector<int>> scores;
	il.init(std::string(argv[1]),'#');
	for(auto &l:il){
		if(l[0]=="scorer") {
			for(int m=1;m<l.size();++m){
				if(nscores.find(l[m])==nscores.end()){
					nscores.insert(std::make_pair(l[m],0));
				}
				nscores.at(l[m])++;
			}
		} else {
			if(scores.find(l[0])==scores.end())
				scores.insert(std::make_pair(l[0],std::vector<int>()));
			for(int m=1;m<l.size();++m)
				scores.at(l[0]).push_back(std::stoi(l[m]));
		}
	}
	for(auto & ns:nscores){
		std::cout << ns.first <<" "<<ns.second<<std::endl;
	}
	for(auto &ss:scores){
		double av=0.0;
		for(auto sc:ss.second) av+=sc;
		av /= (double) ss.second.size();
		std::cout << ss.first<<" " <<av<<std::endl;
	}
}


