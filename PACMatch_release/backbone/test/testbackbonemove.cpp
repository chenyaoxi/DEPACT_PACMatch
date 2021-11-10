/*
 * testbackbonemove.cpp
 *
 *  Created on: 2017年8月1日
 *      Author: hyliu
 */
#include "backbone/backbonemoves.h"
#include "dstl/randomengine.h"
#include <iostream>
#include <fstream>

using namespace NSPproteinrep;

int main(int argc,char **argv){
	std::vector<BackBoneSite> chain;
	readbackbonesites(std::string(argv[1]),chain);
	int seed=std::stoi(std::string(argv[2]));
	NSPdstl::RandomEngine<>::getinstance().reseed(seed);
	BackBoneMoveSelector ms(chain.size());
	BackBoneMoves moves;
	ms.selectmoves(chain,&moves);
	int nloops=moves.nloops();
	for(int i=0; i<nloops;++i) {
		moves.updatechain(&chain,false,i);
		std::string filename="moved"+std::to_string(i)+std::to_string(seed)+".pdb";
		std::ofstream ofs;
		ofs.open(filename.c_str());
		writeSitesToPDB(ofs,chain);
		ofs.close();
	}
}

