/*
 * sdminimizer.cpp
 *
 *  Created on: 2017年11月7日
 *      Author: hyliu
 */
#include "sd/sdrun.h"
#include "backbone/backboneenergy.h"
#include "dstl/randomengine.h"

using namespace NSPproteinrep;
using namespace NSPsd;
double mctemperature(int mcstep, int *minsdsteps,int *maxsdsteps,double *meanvp){
	int cycle=6000;
	mcstep = mcstep%cycle;
	if(mcstep<1000) {
		*minsdsteps=50;
		*maxsdsteps=200;
		*meanvp=0.36;
		return 0.05;
	} else {
		*minsdsteps=5;
		*maxsdsteps=100;
		*meanvp=0.618;
		return 0.005;
	}
}
std::vector<StructRestraint> make_ss_restraints(const std::vector<BackBoneSite> & chain){
	std::vector<double> crd=extractcrd(chain);
	std::vector<std::pair<int,int>> posilen_site;
	std::vector<int> elemid;
	sselements(chain,&elemid,&posilen_site);
	for(auto &pl:posilen_site){
		pl.first *=4;   //every sites has 4 atoms
		pl.second *=4;
	}
	for(auto &c:crd) c*=A2NM;
	std::vector<StructRestraint> results;
	double kres=10000.0;
	for(auto &pl:posilen_site){
		results.push_back(StructRestraint(crd,std::vector<std::pair<int,int>>(1,pl),kres));
	}
	return results;
}
void choosefixatoms(int chainsize,double vportion,std::vector<bool> *fixatoms){
	auto & rng=NSPdstl::RandomEngine<>::getinstance();
	fixatoms->assign(4*chainsize,true);
	int imove=rng.intrng(0,chainsize-1)();
	rng.setrealrng(0,1.0);
	for(int i=0;i<chainsize-3;++i) {
		double p=rng.realrng()();
		if(i==imove || p>vportion){
			for(int idx=4*i; idx<4*i+12;++idx){
				fixatoms->at(idx)=false;
			}
		}
	}
}

int main(int argc,char ** argv) {
	std::vector<BackBoneSite> chain;
	readbackbonesites(std::string(argv[1]), chain);
	std::vector<StructRestraint> ssrestraints=make_ss_restraints(chain); //need "HEC" representation of SS in  chain
	std::vector<std::string> energycontrolines {
		{"PhiPsiWeight=0.0"},
		{"ClashWeight=0.0"},
		{"PBLocalWeight=1.0" },
		{ "PBPackingWeight=0.3" },
		{"RefPBSeq="}
	};
	std::string parasetname1("paraset1");
	adjustenergycontrols(parasetname1,energycontrolines);
	ChainEnergyControl ce1=prepareenergycontrol(&chain,parasetname1);
	BackBoneEnergy energy;
	double eold = energy.totalenergy(chain, &ce1);
	std::vector<double> initcrd=extractcrd(chain);
	for(auto &c:initcrd) c *=A2NM;
	int nmcsteps=std::stoi(std::string(argv[2]));
	int maxsdsteps;
	int minsdsteps;
	double meanvp;
	double mctemp_scale=mctemperature(0,&maxsdsteps,&minsdsteps,&meanvp);
	std::vector<bool> fixatoms;
	choosefixatoms(chain.size(),meanvp,&fixatoms);
	SDRun::SDRunIn sdrunin(initcrd,fixatoms);
	unsigned int seed=31u;
	SDRun sdrun=make_backbone_sdrun(sdrunin,seed);
	for(auto &ssres:ssrestraints) sdrun.ff()->addstructrestraint(ssres);
	int nsdsteps=minsdsteps+maxsdsteps*(0.8*0.8*0.8);
	bool rejected=false;
	double enew;
	double varsd=0.36;
	double mctemp=0.2*(double) chain.size()*mctemp_scale;
	int nstepaccum=0;
	for (int i=0;i<nmcsteps;++i) {
		if( !sdrun.runsteps(nsdsteps)){
			std::cout <<"SHAKE FAILURE"<<std::endl;
			rejected=true;
		} else {
			std::vector<double> crd=sdrun.state().crd;
			for(auto &c:crd) c /=A2NM;
			assigncrd(crd,chain);
			enew=energy.totalenergy(chain,&ce1);
			double de=enew-eold;
			if(de<0.0) rejected=false;
			else if(de <1.e-6){
				rejected=true; //no atoms are moved
			} else{
				rejected=true;
				double paccept=paccept = exp(-de / mctemp);
				if (NSPdstl::RandomEngine<>::getinstance().realrng(0.0,1.0)()< paccept){
					rejected=false;
				}
			}
			std::cout <<"mcstep: " <<i <<" Ene: " <<eold <<" "<<de
					<<" mctemp: " <<mctemp
					<<" nsdsteps: " <<nsdsteps<<std::endl;
		}
		if(rejected ||nstepaccum>10) {
			nstepaccum=0;
			double varp=NSPdstl::RandomEngine<>::getinstance().randomnormal(meanvp,varsd);
			if(varp>0.8) varp=0.8;
			if(varp<0.0) varp=0.0;
			mctemp=(1-varp)*(double)chain.size()*mctemperature(i,&minsdsteps,&maxsdsteps,&meanvp);
			choosefixatoms(chain.size(),varp,&fixatoms);
			sdrun.initstate(SDRun::SDRunIn(initcrd,fixatoms));
//			mctemp=mctemp_scale;
			nsdsteps=minsdsteps+maxsdsteps*(varp*varp*varp);
		}  else {
			initcrd=sdrun.state().crd;
			++nstepaccum;
			eold=enew;
		}

		if((i+1) %100 == 0) {
			for(auto &c:initcrd) c /=A2NM;
			assigncrd(initcrd,chain);
			std::ofstream ofs;
			std::string filename="sdout.pdb";
			ofs.open(filename.c_str());
			writeSitesToPDB(ofs,chain);
			ofs.close();
			for(auto &c:initcrd) c *=A2NM;
		}
	}
}



