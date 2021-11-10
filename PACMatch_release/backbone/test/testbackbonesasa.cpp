/*
 * testbackbonesasa.cpp
 *
 *  Created on: 2016年12月10日
 *      Author: hyliu
 */
#include "backbone/backbonesasa.h"
#include "backbone/segments.h"
#include "backbone/hbonded.h"
#include <iostream>
#include <fstream>
using namespace NSPproteinrep;

int main(int argc, char **argv) {
	std::ifstream ifs;
	ifs.open(argv[1]); // segments files
	std::vector<std::vector<BackBoneSite>> segments;
	readsegments(ifs,segments);
	ifs.close();
	std::vector<std::vector<BackBoneSASA>> sasasegments;
	for( auto & seg:segments){
		sasasegments.push_back(std::vector<BackBoneSASA>());
		std::vector<BackBoneSASA> & sasa=sasasegments.back();
		segmentSASA(seg,sasa);
	}
	for(unsigned int i=0;i<sasasegments.size()-1;++i){
		for(unsigned int j=i+1; j<sasasegments.size();++j) {
			updatesegmentSASA(sasasegments[i],sasasegments[j]);
		}
	}
	std::vector<std::vector<unsigned int>> hbonded;
	hbonded.resize(segments.size(),std::vector<unsigned int>());
	for(unsigned int i=0; i<segments.size();++i){
		hbonded[i].resize(segments[i].size(),0u);
	}
	for(unsigned int i=0;i<segments.size();++i){
		for(unsigned int m=0;m<segments[i].size()-1;++m){
			for(unsigned int n=m+1;n<segments[i].size();++n){
				unsigned int hb=NSPproteinrep::hbonded(segments[i][m],segments[i][n]);
				if(hb==1) {
					hbonded[i][m]+=1;
					hbonded[i][n]+=2;
				}				else if(hb==2){
					hbonded[i][m]+=2;
					hbonded[i][n]+=1;
				}
			}
		}
		if(i==segments.size()-1) continue;
		for(unsigned int j=i+1;j<segments.size();++j) {
				for(unsigned int m=0;m<segments[i].size();++m){
					for(unsigned int n=0;n<segments[j].size();++n){
						unsigned int hb=NSPproteinrep::hbonded(segments[i][m],segments[j][n]);
						if(hb==1) {
							hbonded[i][m]+=1;
							hbonded[j][n]+=2;
						}				else if(hb==2){
							hbonded[i][m]+=2;
							hbonded[j][n]+=1;
						}
					}
				}
		}
	}

	for(unsigned int i=0;i<sasasegments.size();++i){
		for(unsigned int j=0; j<sasasegments[i].size();++j){
			int resid=segments[i][j].resid;
			std::string resname=segments[i][j].resname;
			std::vector<double> exposed;
			sasasegments[i][j].getexposed(&exposed);
			std::cout <<resid<<resname<<" HBOND "<<hbonded[i][j]<< " N: " <<exposed[0] <<std::endl;
			std::cout <<resid<<resname<<" HBOND "<<hbonded[i][j]<<" CA: " <<exposed[1] <<std::endl;
			std::cout <<resid<<resname<<" HBOND "<<hbonded[i][j]<<" C: " <<exposed[2]<<std::endl;
			std::cout <<resid<<resname<<" HBOND "<<hbonded[i][j]<<" O: " <<exposed[3]<<std::endl;
		}
	}
}


