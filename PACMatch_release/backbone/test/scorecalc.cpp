/*
 * scorecalc.cpp
 *
 *  Created on: 2017年1月3日
 *      Author: hyliu
 */
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include "dataio/inputlines.h"
int main(int argc, char**argv){
        NSPdataio::InputLines lines;
        lines.init(std::string(argv[1]),'#');
        NSPdataio::InputLines lines2;
        lines2.init(std::string(argv[2]),'#');
        double w1=std::stod(std::string(argv[3]));
        double diff=0.0;
        double gap=0.0;
        double np=0;
        double nn=0;
        unsigned int idx2=0;
        for(auto & l:lines) {
        		std::vector<std::string> &l2=lines2[idx2++];
        		assert(l[0] == l2[0]);
                double score_native;
                score_native=w1*std::stod(l[1]) + std::stod(l2[1]);
                std::vector<double> scores(10,0.0);
                double min=1000000000.0;
                for(unsigned int i=0;i<10;++i){
                        scores[i]=w1*std::stod(l[i+2])+std::stod(l2[i+2]);
                        if(scores[i]<min)min=scores[i];
                }
                double sav=0.0;
                double sfluc=0.0;
                for(auto s:scores) {
                	sav +=s;
                	sfluc+=s*s;
                }
                sav /=10.0;
                sfluc= sqrt(sfluc/10.0-sav*sav);
                gap +=(score_native-sav)/sfluc;
   //             std::cout <<min <<"\t" <<score_native <<std::endl;
                if(min > score_native) ++np;
                else ++nn;
        }
        std::cout <<"score_gap " <<gap/(double) (lines.size()) <<std::endl;
        std::cout <<"positive_ratio: "<<np/(np+nn) <<std::endl;
}

