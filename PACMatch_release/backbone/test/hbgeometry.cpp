/*
 * hbgeometry.cpp
 *
 *  Created on: 2016年5月13日
 *      Author: hyliu
 */



#include "backbone/backbonesite.h"
#include "geometry/relativeposition.h"
#include "backbone/hbonded.h"
#include <iostream>
#include <fstream>
using namespace NSPgeometry;
using namespace NSPproteinrep;
//#include "geometry/histogram1d.h"
int main(int argc, char **argv) {
	if (argc < 3) {
		std::cout << "usage: " << argv[0] << " sitefile pairfile\n";
		exit(0);
	}
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
	std::ifstream ifs;
	ifs.open(argv[2]);
	long npairs=std::stol(std::string(argv[3]));
	unsigned int nhbonded=0;
	for(long i=0; i<npairs;++i) {
		long i1;
		long i2;
		ifs >>i1;
		ifs>>i2;
		RelativePosition rp;
		rp.read(ifs);
		BackBoneSite &s1=sites[i1];
		BackBoneSite &s2=sites[i2];
		HBondGeometry hbg;
		if(hbonded(s1,s2,&hbg)){
			nhbonded++;
			std::cout <<hbg.rha <<"\t" <<hbg.ahab <<"\t" <<hbg.adha<<"\t"
					<<hbg.rda <<"\t"<< hbg.adab<< std::endl;
		}
	}
	std::cout <<nhbonded <<std::endl;
}

