/*
 * checksitesfile.cpp
 *
 *  Created on: 2017年5月24日
 *      Author: hyliu
 */

#include "backbone/mainchain.h"
using namespace NSPproteinrep;

#include "backbone/mainchain.h"

int main(int argc, char **argv) {
	MainChain chain;
	chain.read(std::string(argv[1]));
	for(int i=0; i<chain.size();++i){
		double phi,psi,omiga;
		double rphi{0},rpsi{0},romiga{0};
		phi=chain[i].phi();
		psi=chain[i].psi();
		omiga=chain[i].omiga();
		if(i>0) {
			rphi=chain[i].phi(chain[i-1]);
		}
		if(i<chain.size()-1){
			rpsi=chain[i].psi(chain[i+1]);
			romiga=chain[i].omiga(chain[i+1]);
		}
		std::cout <<i<<"\t("<<phi<<","<<psi<<","<<omiga<<")\t(";
		std::cout <<rphi<<","<<rpsi<<","<<romiga<<")"<<std::endl;
	}
}
