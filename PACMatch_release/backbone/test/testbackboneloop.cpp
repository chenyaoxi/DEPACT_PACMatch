/*
 * testbackboneloop.cpp
 *
 *  Created on: 2016年12月12日
 *      Author: hyliu
 */
#include "backbone/backboneloop.h"
#include "backbone/loopinteractions.h"
#include <iostream>
#include <fstream>
using namespace NSPproteinrep;

/**
 * This test program rebuilt loops for a series of pdb structures contained in the backbone site file
 * specified in argv[1]. Locations of the loops are provided in the file specified by
 *  argv[2]. For each loop, ten structures will be built. The built loops are written to
 *  "generatedloops.dat".
 */
int main(int argc, char **argv) {
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
	std::vector<std::vector<long>> loops;
	std::vector<double> refscores;
	std::ifstream ifs;
	ifs.open(argv[2]); //loop locations and reference scores
	while (ifs.good()) {
		long l0, l1, l2, l3;  //l0 start of the chain, l3 end of the chain
							  //l1 start of the loop, l2 end of the loop
		double rs;
		ifs >> l0 >> l1 >> l2 >> l3 >>rs;
		if (ifs.good()) {
			loops.push_back(std::vector<long>( { l0, l1, l2, l3 }));
			refscores.push_back(rs);
		}
	}
	unsigned int nloop = 0;
	IdealGeometries & idg = IdealGeometries::getGlobalInstance(
			"idealgeometries.dat");
	std::ofstream ofs;
	ofs.open("generatedloops.dat");
	LoopInteractions::switches.tetrasef_on=false;
	LoopInteractions::switches.torsionmotif_on=false;
	for (auto & lp : loops) {
		if(nloop <=804) {
			++nloop;
			continue;
		}
		std::vector<std::vector<BackBoneSite>> segments;
		std::vector<BackBoneSite> loop;
		segments.push_back(std::vector<BackBoneSite>());
		segments[0].resize(lp[1] - lp[0], BackBoneSite()); //contains partial chain ahead of the loop
		std::copy(sites.begin() + lp[0], sites.begin() + lp[1],
				segments[0].begin());
		segments.push_back(std::vector<BackBoneSite>());
		segments[1].resize(lp[3] - lp[2], BackBoneSite()); //partial chain after the loop
		std::copy(sites.begin() + lp[2], sites.begin() + lp[3],
				segments[1].begin());
		loop.resize(lp[2] - lp[1], BackBoneSite());
		std::copy(sites.begin() + lp[1], sites.begin() + lp[2], loop.begin());
		BackBoneLoop rloop(&segments, 0u, 1u, loop.size(), loop); //loop is used here to
		//inform the constructor of gly and pro locations
		rloop.setsaveacceptable(refscores[nloop]*0.99);
		unsigned int ncopies = 10;
		unsigned int ncandidates = 10;
		std::vector<std::vector<BackBoneSite>> newloops;
		newloops.resize(ncopies + 1, std::vector<BackBoneSite>());
		//newloops is allocated one-loop larger, the first position could be
		//used to hold the native loop (not used here).
		std::vector<double> scores(ncopies,0);
		ncopies = rloop.getLoops(ncopies, ncandidates, newloops.begin() + 1,&scores);
		//The generated loops are stored in newloops, from index 1.
		ofs << nloop << " " << ncopies << std::endl;
		for(int i=0; i<ncopies;++i) {
			std::cout <<scores[i]<<"\t";
		}
		std::cout <<std::endl;
		for (unsigned int m = 0; m < ncopies; m++) {
			for (auto &bs : newloops[m + 1])
				ofs << bs.toString();
		}
		ofs.flush();
		/*
		 newloops[0].resize(loop.size()+2, BackBoneSite());
		 std::copy(sites.begin() + lp[1] - 1, sites.begin() + lp[2] + 1,
		 newloops[0].begin());
		 segments.clear();
		 segments.push_back(std::vector<BackBoneSite>());
		 segments[0].resize(lp[1] - lp[0] - 1, BackBoneSite());
		 std::copy(sites.begin() + lp[0], sites.begin() + lp[1] - 1,
		 segments[0].begin());
		 segments.push_back(std::vector<BackBoneSite>());
		 segments[1].resize(lp[3] - lp[2] - 1, BackBoneSite());
		 std::copy(sites.begin() + lp[2] + 1, sites.begin() + lp[3],
		 segments[1].begin());

		 for(unsigned int i=0;i<ncopies+1;++i) {
		 std::vector<BackBoneSite> chain;
		 chain.resize(segments[0].size()+newloops[0].size()+segments[1].size(), BackBoneSite());
		 std::copy(segments[0].begin(),segments[0].end(),chain.begin());
		 std::copy(newloops[i].begin(),newloops[i].end(),chain.begin()+segments[0].size());
		 std::copy(segments[1].begin(),segments[1].end(),chain.begin()+segments[0].size()+newloops[0].size());
		 std::string filename="loop"+std::to_string(nloop)+"_"+std::to_string(i)+".pdb";
		 std::ofstream ofs;
		 ofs.open(filename.c_str());
		 writeSitesToPDB(ofs,chain);
		 ofs.close();
		 }*/
		++nloop;
	}
	/*
	 std::vector<std::vector<BackBoneSASA>> segmentsasa;
	 std::vector<std::vector<unsigned int>> segmenthbonded;
	 for (auto & seg : segments) {
	 segmentsasa.push_back(std::vector<BackBoneSASA>());
	 segmentSASA(seg, segmentsasa.back());
	 segmenthbonded.push_back(std::vector<unsigned int>());
	 NSPproteinrep::segmenthbonded(seg,segmenthbonded.back());
	 }




	 updatesegmentSASA(segmentsasa[0], segmentsasa[1]);
	 updatehbonded(segments[0],segmenthbonded[0],segments[1],segmenthbonded[1]);
	 for (unsigned int i=0; i<ncopies + 1; ++i) {
	 std::vector<std::vector<BackBoneSASA>> sasa_f=segmentsasa;
	 std::vector<std::vector<unsigned int>> hbonded_f=segmenthbonded;
	 std::vector<BackBoneSASA> loopsasa;
	 std::vector<unsigned int> loophbonded;
	 segmentSASA(newloops[i],loopsasa);
	 NSPproteinrep::segmenthbonded(newloops[i],loophbonded);
	 updatesegmentSASA(sasa_f[0], loopsasa);
	 updatesegmentSASA(sasa_f[1], loopsasa);
	 updatehbonded(segments[0],hbonded_f[0],newloops[i],loophbonded);
	 updatehbonded(segments[1],hbonded_f[1],newloops[i],loophbonded);
	 for(unsigned int s=0; s<2; ++s){
	 for(unsigned int p=0; p<segments[i].size;++p) {
	 std::vector<double> exposed;
	 std::vector<double> exposed_c;
	 sasa_f[s][p].getexposed(&exposed);
	 segmentsasa[s][p].getexposed(&exposed_c);
	 if(exposed[0] != exposed_c[0]) {
	 std::cout <<i<<" "<<"N_expose_change: " <<segmenthbonded[s][p]%10 <<" " <<hbonded_f[s][p]%10 <<" "
	 <<exposed_c[0] <<" "<< exposed[0] <<std::endl;
	 }
	 if(exposed[3] != exposed_c[3]) {
	 std::cout <<i<<" "<<"C_expose_change: " <<segmenthbonded[s][p]/10u <<" " <<hbonded_f[s][p]/10u <<" "
	 <<exposed_c[3] <<" "<< exposed[3] <<std::endl;
	 }
	 }
	 }
	 for(unsigned int p=0;p<newloops[i].size();++p) {
	 std::vector<double>exposed;
	 loopsasa[p].getexposed(&exposed);
	 unsigned int nhbonded=loophbonded[i]%10;
	 unsigned int ohbonded=loophbonded[i]/10u;
	 if(newloops[i][p].resname != "PRO")
	 std::cout <<i<<" "<<"N: " << nhbonded <<" " <<exposed[0] <<std::endl;
	 std::cout <<i<<" "<<"O: " << ohbonded <<" " <<exposed [3]  <<std::endl;
	 }

	 }
	 }
	 */
}

