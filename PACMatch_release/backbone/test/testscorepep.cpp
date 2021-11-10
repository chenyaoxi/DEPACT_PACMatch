/*
 * testscorepep.cpp
 *
 *  Created on: 2017年6月30日
 *      Author: hyliu
 */
#include "backbone/scorepep.h"
#include "pdbstatistics/proteinblock.h"
#include "pdbstatistics/phipsidistr.h"
using namespace NSPproteinrep;
using namespace NSPgeometry;

int main(int argc,char **argv){
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
	std::vector<std::vector<double>> data;
	std::vector<int> positions;
//	std::vector<PepConformer> conformers;t
	int conftype=std::stoi(std::string(argv[2]));
	ScorePep scorepep(conftype);
	int peptype=std::stoi(std::string(argv[3]));
	int peplength=5;
	scorepep.buildtree(sites,peplength);
//	scorepep.buildreftree(5*scorepep.templatesize());
	scorepep.buildrefphipsitrees(50000);
	std::ofstream ofs_n;
	std::ofstream ofs_r;
	std::string filename1;
	std::string filename2;
	if(peptype==ScorePep::COIL){
		filename1="ncoilc_ene.dat";
		filename2="rcoilc_ene.dat";
	} else if(peptype==ScorePep::HELIX){
		filename1="nhelixc_ene.dat";
		filename2="rhelixc_ene.dat";
	} else if(peptype==ScorePep::STRAND){
		filename1="nstrandc_ene.dat";
		filename2="rstrandc_ene.dat";
	} else {
		filename1="nativec_ene.dat";
		filename2="randomc_ene.dat";
	}
	ofs_n.open(filename1.c_str());
	ofs_r.open(filename2.c_str());
	int nsample=0;
	double sx,sr_ana;
	const NSPpdbstatistics::PhiPsiDistr *distr=&(NSPpdbstatistics::PhiPsiDistr::mixcoildistr());
	for(long i=2;i<sites.size()-2; ++i) {
		if(nsample > 110000) break;
		std::vector<double> conf=scorepep.extractconf(sites,i,peplength);
		if(conf.empty()) continue;
		char pbtype=NSPpdbstatistics::ProteinBlock::pbtype(conf);
		if(peptype ==ScorePep::HELIX){
			if(pbtype != 'm' || sites[i].sscodechar() != 'H') continue;
		} else if(peptype==ScorePep::STRAND){
			if(pbtype != 'd' || sites[i].sscodechar() != 'E') continue;
		} else if(peptype==ScorePep::COIL){
			if(sites[i].sscodechar() != 'C' &&(pbtype=='m' || pbtype=='d')) continue;
		}
		double d2m,d20;
		double s= -log(scorepep.score(conf,&sx,&sr_ana,&d2m,&d20));
		double ephipsi=distr->psi_statisticalenergy(conf[0])+
				distr->statisticalenergy(conf[1],conf[2])+distr->statisticalenergy(conf[3],conf[4])+
				distr->statisticalenergy(conf[5],conf[6])+distr->phi_statisticalenergy(conf[7]);
//		double s=scorepep.neighborsum(conf,scorepep.tree());
//		std::cout <<s<<' '<<log(s)<<":";
//		std::cout <<std::endl;
		for (auto d:conf) ofs_n <<d <<' ';
		ofs_n <<s <<" "<<s+ephipsi<<"  "<<sx<<"  "<<sr_ana<<"  "<<d2m<<" "<<d20<<" "<<pbtype<<std::endl;
		std::vector<double> conf_r=scorepep.sampleconf(peplength,peptype);
		pbtype=NSPpdbstatistics::ProteinBlock::pbtype(conf_r);
		double ephipsi_r=distr->psi_statisticalenergy(conf_r[0])+
					distr->statisticalenergy(conf_r[1],conf_r[2])+distr->statisticalenergy(conf_r[3],conf_r[4])+
					distr->statisticalenergy(conf_r[5],conf_r[6])+distr->phi_statisticalenergy(conf_r[7]);
		double sr=-log(scorepep.score(conf_r,&sx,&sr_ana,&d2m,&d20));
//		double sr=scorepep.neighborsum(conf_r,scorepep.tree());
		for (auto d:conf_r) ofs_r <<d <<' ';
		ofs_r <<sr<<" "<<sr+ephipsi_r<<"  "<<sx<<"  "<<sr_ana<<" "<<d2m<<" "<<d20<<" "<<pbtype<<std::endl;
		++nsample;
	}
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
