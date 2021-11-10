/*
 * gaminimizer.cpp
 *
 *  Created on: 2017年11月13日
 *      Author: hyliu
 */
#include "sd/sdrun.h"
#include "backbone/backboneenergy.h"
#include "dstl/randomengine.h"
#include "dstl/ga.h"
#include "geometry/quatfit.h"
#include <thread>
#include <mutex>
using namespace NSPproteinrep;
using namespace NSPsd;
std::vector<StructRestraint> make_ss_restraints(
		const std::vector<BackBoneSite> & chain) {
	std::vector<double> crd = extractcrd(chain);
	std::vector<std::pair<int, int>> posilen_site;
	std::vector<int> elemid;
	sselements(chain, &elemid, &posilen_site);
	for (auto &pl : posilen_site) {
		pl.first *= 4;   //every sites has 4 atoms
		pl.second *= 4;
	}
	for (auto &c : crd)
		c *= A2NM;
	std::vector<StructRestraint> results;
	double kres = 10000.0;
	for (auto &pl : posilen_site) {
		results.push_back(
				StructRestraint(crd, std::vector<std::pair<int, int>>(1, pl),
						kres));
	}
	return results;
}
/*void choosefixatoms(int chainsize,double vportion,std::vector<bool> *fixatoms){
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
 }*/
struct SDThreadArg {
	int thread_id;
	SDRun *thread_sdrun;
	ChainEnergyControl *ce;
	BackBoneEnergy *energy;
	std::vector<BackBoneSite> *thread_chain;
	std::vector<double> *initcrd;
	int nsdsteps;
	std::mutex *gamutex;
	NSPdstl::GAPopulation<std::vector<double>> *gapopulation;

};
/**
 * compare a set of coordinate crd2 to a reference set crd1
 *  and generate a score contributing to fitness for GA ranking of crd2.
 *  If crd1 and crd2 point to the same data or rmsd2 between crd2 and crd1 is larger than rmsdcut2, a minscore(best fitness)
 *  is returned.
 *  Otherwise, the score is a large negative number times rmsd2,so that the final fitness will
 *  be mainly come from rmsd2.
 *
 */
class RmsdScoreForGA {
public:
	double operator()(std::shared_ptr<std::vector<double>> crd1,
			std::shared_ptr<std::vector<double>> crd2) const {
		const static double rmsdcut2 = 0.0025;
		const static double scale = -1.e8;
		static double scoremin = scale * rmsdcut2;
		if (crd1 == crd2)
			return scoremin;
		double rmsd2 = NSPgeometry::QuatFit().setup(*crd1, *crd2);
		if (rmsd2 < rmsdcut2)
			return scale * rmsd2; //larger rmsd,lower score, better fitness
		return scoremin;
	}
};

void threadsdrunner(SDThreadArg *arg) {
	SDThreadArg *targ = (SDThreadArg *) arg;
	targ->thread_sdrun->initstate(SDRun::SDRunIn(*(targ->initcrd)));
	if (!targ->thread_sdrun->runsteps(targ->nsdsteps)) {
		std::cout << "SHAKE FAILURE!!!!" << std::endl;
		return;
	}
	std::shared_ptr<std::vector<double>> newcrd = std::shared_ptr<
			std::vector<double>>(
			new std::vector<double>(targ->thread_sdrun->state().crd));
	for (auto &c : *newcrd)
		c /= A2NM;
	assigncrd(*newcrd, *(targ->thread_chain));
//	std::cout <<"threaid "<<targ->thread_id <<(*(targ->thread_chain))[0].ncrd().x_<<" "
//			 <<	(*(targ->thread_chain)).back().ncrd().x_<<std::endl;
	double enew = targ->energy->totalenergy(*(targ->thread_chain), targ->ce);
//	std::cout <<"threaid "<<targ->thread_id <<"New ene "<< enew	<<std::endl;
	for (auto &c : *newcrd)
		c *= A2NM;
	std::lock_guard < std::mutex > lock(*(targ->gamutex));
	targ->gapopulation->addindividual(newcrd, enew);
	return;
}
int main(int argc, char ** argv) {
	std::vector<BackBoneSite> chain;
	readbackbonesites(std::string(argv[1]), chain);
	std::vector<StructRestraint> ssrestraints = make_ss_restraints(chain); //need "HEC" representation of SS in  chain
	std::vector<std::string> energycontrolines { { "PhiPsiWeight=0.0" }, {
			"ClashWeight=0.0" }, { "PBLocalWeight=1.0" }, {
			"PBPackingWeight=0.3" }, { "RefPBSeq=" },{"IgnoreResName=1"} };
	int maxgenerations = 2000;
	int maxstored = 1000;
	int nchildren = 3;
	int long_sdsteps = 599;
	int short_sdsteps=5;
	int popsize = 14;
	int subsize = 3;
	int nthread = 2;
	std::vector<unsigned int> seeds;
	for (int i = 0; i < nthread + 1; ++i)
		seeds.push_back(((i + 13) * 311) % 137 + 7);
	std::string parasetname1("paraset1");
	adjustenergycontrols(parasetname1, energycontrolines);
	ChainEnergyControl ce1 = prepareenergycontrol(&chain, parasetname1);
	BackBoneEnergy energy(ce1);  //ignoreresname in energy_control will be passed
	double eold = energy.totalenergy(chain, &ce1);
	std::shared_ptr<std::vector<double>> initcrd = std::shared_ptr<
			std::vector<double>>(new std::vector<double>(extractcrd(chain)));
	for (auto &c : *initcrd)
		c *= A2NM;
	SDRun::SDRunIn sdrunin(*initcrd);
	SDRun sdrun = make_backbone_sdrun(sdrunin, seeds.back());
	for (auto &ssres : ssrestraints)
		sdrun.ff()->addstructrestraint(ssres);
	NSPdstl::GAPopulation<std::vector<double>> gapopulation(popsize, subsize);
	gapopulation.addindividual(initcrd, eold);
	RmsdScoreForGA rmsdscorer;
	assert(seeds.size() >= nthread);
	std::vector<BackBoneEnergy> threadenergy(nthread, energy);
	std::vector<SDRun> threadsdruns(nthread, sdrun);
	std::vector<ChainEnergyControl> threadce(nthread, ce1);
	std::vector<std::vector<BackBoneSite>> threadchain(nthread, chain);
	std::vector<SDThreadArg> threadargs(nthread);
	std::mutex gamutex;
	for (int i = 0; i < nthread; ++i) {
		threadargs[i].thread_id = i;
		threadargs[i].thread_sdrun = &(threadsdruns[i]);
		threadargs[i].thread_sdrun->initrandomengine(seeds[i]);
		threadargs[i].ce = &(threadce[i]);
		threadargs[i].energy = &(threadenergy[i]);
		threadargs[i].thread_chain = &(threadchain[i]);
		threadargs[i].gamutex = &gamutex;
		threadargs[i].gapopulation = &gapopulation;
//		threadargs[i].nsdsteps = nsdsteps;
	}
	std::vector < std::shared_ptr
			< std::thread >> startedthreads(nthread, nullptr);
	int ngenafterstore = 0;
	for (int gen = 0; gen < maxgenerations; ++gen) {
		auto & ivs = gapopulation.individuals();
		int nivs = ivs.size();
		int nruns = 0;
		int nrunning = 0;
		for (int m = 0; m < nivs; ++m) {
			int nc=1;
			if(m==0||m<gapopulation.ncenter()*subsize){
				nc=nchildren;
			}
			for (int c = 0; c < nc; ++c) {
				int myid = nruns % nthread;
				if (nruns >= nthread) {
					startedthreads[myid]->join();
					--nrunning;
				}
				if(m==0||m<gapopulation.ncenter()*subsize){
					threadargs[myid].nsdsteps=short_sdsteps;
				} else{
					threadargs[myid].nsdsteps=long_sdsteps;
				}
				threadargs[myid].initcrd = ivs.at(m).get();
				startedthreads[myid] = std::shared_ptr < std::thread
						> (new std::thread(threadsdrunner, &threadargs[myid]));
//				startedthreads[myid]->join();
				++nrunning;
				++nruns;
			}
		}
		for (int i = 0; i < nrunning; ++i) {
			startedthreads[i]->join();
		}
		int ncenters = gapopulation.selection_distance(rmsdscorer);
		//print energies
		const std::vector<double> &fitnesses = gapopulation.fitness();
		std::cout << "Generation: " << gen <<"; Population_size: "
				<< fitnesses.size() << "; Number_Pop_Centers: " << ncenters
				<< std::endl;
		for (auto f : fitnesses)
			std::cout << "\t" << f;
		std::cout << std::endl;
		if (ncenters * subsize >= popsize) {
			++ngenafterstore;
			if (ngenafterstore >= 5) {
				ngenafterstore = 0;
				gapopulation.StoreNShrink(rmsdscorer);
				auto &stored = gapopulation.visited();
				//write coordinates as pdb files
				for (int m = 0; m < stored.size(); ++m) {
					std::vector<double> newcrd = *(stored.at(m).first);
					for (auto &c : newcrd)
						c /= A2NM;
					assigncrd(newcrd, chain);
					std::ofstream ofs;
					std::string filename = "garesult" + std::to_string(m)
							+ ".pdb";
					ofs.open(filename.c_str());
					ofs << "Energy: " << stored.at(m).second << std::endl;
					writeSitesToPDB(ofs, chain);
					ofs.close();
				}
				if (stored.size() >= maxstored)
					break;
			}
		}
	}
}
