/*
 * decomposesiteenergies.cpp
 *
 *  Created on: 2017年5月26日
 *      Author: hyliu
 */


#include "backbone/mainchain.h"
#include "backbone/energyfunctions.h"

using namespace NSPproteinrep;
int main(int argc, char **argv) {
	EnergyTerms erank;
	erank.addtorsionvecene();
	erank.addtetrasefene(0.17);
	MainChain chain;
	chain.read(std::string(argv[1]));
	decomposeenergy(chain,erank,std::cout);
}
