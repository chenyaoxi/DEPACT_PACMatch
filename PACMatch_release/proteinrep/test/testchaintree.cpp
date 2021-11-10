/*
 * testchaintree.cpp
 *
 *  Created on: 2016年11月23日
 *      Author: hyliu
 */

#include "proteinrep/chaintree.h"
#include "proteinrep/pdbrecord.h"
using namespace NSPproteinrep;
int main() {
	std::shared_ptr<ChainTreeTopology> topo=makeMainChainHeavyTopo(20);  //A 6 unit backbone topology
	IdealGeometries & idg=IdealGeometries::getGlobalInstance("idealgeometries.dat"); // read in ideal internal coordinates
	ChainTreeCrd crd(topo.get());
	crd.initWithIdealIntCrd(); // set ideal internal coordinates (they should include initial values  even for rotable torsions)
	crd.initCrdPosi0();   //setup coordinates of first few atoms that can not be determined for internal coordinates
	const double deg=3.14159265358979323846/180.0;
	for(int p=0; p<20;++p) {
		crd.setBackBoneTorsionAt(p,ChainTreeCrd::BackBoneTorsions(-40*deg,-60*deg)); //helical
	}
	crd.writeIntCrd(std::cout);
	crd.calcXYZ();  //calculate coordinates;
	crd.setBackBoneTorsionAt(10,ChainTreeCrd::BackBoneTorsions(40*deg,120*deg));
	crd.calcXYZ();  //calculate coordinates;
	crd.writePDB(std::cout);
}


