/*
 * testmainchain.cpp
 *
 *  Created on: 2017年4月22日
 *      Author: hyliu
 */
#include "backbone/mainchain.h"
#include "backbone/energyfunctions.h"
#include "backbone/closingloop.h"
#include "dstl/randomengine.h"
#include <memory>
using namespace NSPproteinrep;

int main(int argc, char **argv){
	NSPdstl::RandomEngine<>::getinstance().init(55u);
	EnergyTerms eterms;
	eterms.addtorsionene();
	eterms.addtetrasefene();
//	eterms.addtorsionvecene();
//	eterms.addstericclash();
	MainChain mc;

	mc.generaterandom(20);

	std::map<std::string,double> energies;
	calctotalenergy(mc,eterms,energies,1.e20);
	for(auto & e:energies) {
		std::cout << e.first <<"\t"<<e.second <<std::endl;
	}
	std::cout <<std::endl;

	mc.masksegments(std::vector<std::pair<int,int>>({{5,3}}));
	energies.clear();
	calctotalenergy(mc,eterms,energies,1.e20);
	for(auto & e:energies) {
		std::cout << e.first <<"\t"<<e.second <<std::endl;
	}
	std::cout <<std::endl;

	energies.clear();
	std::vector<BackBoneSite> newloop;
	mc.masksegments(std::vector<std::pair<int,int>>({{1,11}}));
	newloop.resize(3);
	std::copy(mc.begin()+5, mc.begin()+8,newloop.begin());
	calcloopenergy(mc,5,3,newloop,eterms,
			energies,1.e20);
	for(auto & e:energies) {
		std::cout << e.first <<"\t"<<e.second <<std::endl;
	}
	std::cout <<std::endl;

//	writeSitesToPDB(std::cout,mc);
//	mc.setrigidbetween(7,11);
//	mc.setrigidbetween(16,20);
//	mc.insertgap(16,3);

/*	ClosingLoop cl(mc,7,10,ClosingLoop::MUTATESEQ_CONF);
	for (int m=0; m< 15;++m) {
		int nsol=cl.solve();
		for (int i=0; i<nsol; ++i){
			std::vector<BackBoneSite> sol;
			cl.getsitessolution(i,&sol);
			std::ofstream ofs;
			std::string filename="solution" +std::to_string(m)+"_"+std::to_string(i)+".pdb";
			ofs.open(filename.c_str());
			writeSitesToPDB(ofs,sol);
			ofs.close();
		}
	}
	*/
}
