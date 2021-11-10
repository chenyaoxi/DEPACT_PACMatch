/*
 * avcoursescore.cpp
 *
 *  Created on: 2017年1月4日
 *      Author: hyliu
 */

#include <iostream>
#include <vector>
#include <map>
#include <cassert>
#include "dataio/inputlines.h"
int main(int argc, char**argv){
        NSPdataio::InputLines lines;
        lines.init(std::string(argv[1]),'#');
        std::map<std::string,std::vector<double>> sscores;
        unsigned int idx2=0;
        for(auto & l:lines) {
        		std::vector<double> scores;
        		for(int i=1; i<l.size();++i) {
        			double s=std::stod(l[i]);
        			assert (s>50 && s<=101);
        			scores.push_back(s);
        		}
        		sscores.insert(std::make_pair(l[0],scores));
        }
        for( auto &l:lines){
        	std::vector<double> &ss=sscores.at(l[0]);
        	double sav=0.0;
        	for( auto s:ss) {
        		sav +=s;
        	}
        	std::cout << l[0] <<"\t"<<sav/(double) (ss.size())<<std::endl;
        }
}



