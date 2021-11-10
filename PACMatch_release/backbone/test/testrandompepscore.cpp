/*
 * testrandompepscore.cpp
 *
 *  Created on: 2017年7月12日
 *      Author: hyliu
 */

#include "backbone/scorepep.h"
#include "pdbstatistics/proteinblock.h"
using namespace NSPproteinrep;
using namespace NSPgeometry;

int main(int argc,char **argv){
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
	int ngly=0;
	int npro=0;
	ScorePep scorepep;
	int peplength=5;
	long posi=0;
	long nres=0;
	for (auto &s:sites) {
		std::vector<double> conf=scorepep.extractconf(sites,posi,peplength);
		++posi;
	        if(conf.empty()) continue;
		char pbtype=NSPpdbstatistics::ProteinBlock::pbtype(conf);
		if(pbtype =='d' || pbtype =='m') continue;	
		++nres;
		if(s.resname=="GLY") ++ngly;
		if(s.resname=="PRO") ++npro;
	}
	std::cout <<"PGLY: " << (double) ngly/(double) nres;
	std::cout <<"PPRO: " << (double) npro/(double) nres;
//	std::vector<std::vector<double>> data;
//	std::vector<int> positions;
//	std::vector<PepConformer> conformers;
/*	ScorePep scorepep;
	int peplength=5;
//	scorepep.buildtree(sites,peplength);
	scorepep.buildreftree(sites.size());
	ScorePep scorepep2;
	scorepep2.buildreftree(sites.size());
	std::ofstream ofs_n;
	std::ofstream ofs_r;
	int nsample=1000;
	for(long i=0;i<nsample; ++i) {
		std::vector<double> conf_r=ScorePep::sampleconf(peplength);
		double sr=scorepep.neighborsum(conf_r,scorepep.reftree());
		double sr2=scorepep2.neighborsum(conf_r,scorepep2.reftree());
		double asr=scorepep.refprobability(conf_r);
		std::cout<<sr<<" "<<sr2<<" "<< asr <<" " <<sr/asr<<" :";
		for(auto c:conf_r) std::cout <<" "<<c;
		std::cout<<std::endl;
		++nsample;
	} */
/*	std::cout <<"random conformers: " <<std::endl;
	for(long i=0;i<nsample; ++i) {
		std::vector<double> conf=ScorePep::sampleconf(peplength);
		double s= scorepep.score(conf);
		std::cout <<s<<' '<<log(s)<<":";
		std::cout <<std::endl;
		for (auto d:conf) ofs <<d<<' ';
		ofs <<log(s) <<std::endl;
	}*/
}



