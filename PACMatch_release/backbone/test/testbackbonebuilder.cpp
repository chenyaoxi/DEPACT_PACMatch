/*
 * testbackbonebuilder.cpp
 *
 *  Created on: 2018年1月9日
 *      Author: hyliu
 */


#include "backbone/backbonebuilder.h"
#include "dstl/randomengine.h"
#include "fullsite/fullsite.h"
using namespace NSPproteinrep;
using namespace NSPgeometry;
int main(int argc,char **argv){
	if(argc>1){
		int seed=std::stoi(std::string(argv[1]));
		NSPdstl::RandomEngine<>::getinstance().reseed(seed);
	}
	std::vector<std::vector<FullSite>> chains=readfullsitesfrompdb("pdb3l32.ent");
	std::vector<BackBoneSite> bss1=backbone(chains[0]);
	std::vector<BackBoneSite> bss2=backbone(chains[1]);
	std::vector<std::shared_ptr<std::vector<BackBoneSite>>> l12;
	int ntry1=0;
	do { l12=BackBoneBuilder::buildlinkers(4,bss1.back(),bss2[0],
		std::vector<std::pair<int,int>> (),
		std::vector<std::pair<int,int>>(),std::set<int>());
		++ntry1;
		std::cout <<"Try linking " << ntry1<<std::endl;
	}	while(l12.empty() && ntry1<10);
	std::vector<std::vector<FullSite>> results(1);
	std::vector<FullSite> &result=results[0];
	for(auto &s:chains[0]) result.push_back(s);
	for (BackBoneSite s:*(l12[0])) {
		s.resname="GLY";
		s.pdbid='A';
		s.resseq=result.size()+1;
		FullSite fs=make_fullsite(s);
		result.push_back(fs);
	}
	for(auto &s:chains[1]) result.push_back(s);
	std::ofstream os;
	std::string fname="3l32linked.pdb";
	os.open(fname.c_str());
	writetopdb(results,os);
	os.close();
	exit(0);
	std::vector<XYZ> strandstarts;
	for(int i=0;i<10;++i) {
		strandstarts.push_back(XYZ(5*i,0.0,0.0));
	}

	XYZ direction(0,1.0,0);
/*	std::vector<std::pair<int,int>> helixregions;
	helixregions.push_back(std::make_pair(0,20));
	helixregions.push_back(std::make_pair(25,20));
	helixregions.push_back(std::make_pair(50,15));
	helixregions.push_back(std::make_pair(70,15));
	BackBoneSite bsn;
	genbackbonesite(nullptr,false,-120.0,120.0,&bsn);
	bsn.chainid='A';
	bsn.resid=-1;
	bsn.resseq=-1;
	std::vector<BackBoneSite> chain=BackBoneBuilder::buildforwardbackbone(85,
			bsn,helixregions,std::vector<std::pair<int,int>> (), std::set<int>());*/
	std::vector<std::vector<BackBoneSite>> strands(4);
	strands[0]=BackBoneBuilder::buildstrandat(5,strandstarts[0],direction,true);
	strands[1]=BackBoneBuilder::buildstrandat(5,strandstarts[1],direction,false);
	strands[2]=BackBoneBuilder::buildstrandat(5,strandstarts[2],direction,true);
	strands[3]=BackBoneBuilder::buildstrandat(5,strandstarts[3],direction,false);
	std::vector<BackBoneSite> nterm=BackBoneBuilder::buildbackwardbackbone(2,
			strands[2][0],std::vector<std::pair<int,int>>(),
			std::vector<std::pair<int,int>>(),std::set<int>());
	std::vector<BackBoneSite> cterm=BackBoneBuilder::buildforwardbackbone(2,
			strands[3].back(),std::vector<std::pair<int,int>>(),
			std::vector<std::pair<int,int>>(),std::set<int>());
	std::vector<XYZ> uhelixstarts;
	std::vector<XYZ> dhelixstarts;
	for(int i=0;i<5;++i){
		uhelixstarts.push_back(XYZ(9*i,0.0,4.0));
	}
	for(int i=0;i<5;++i){
		dhelixstarts.push_back(XYZ(9*i,0.0,-4.0));
	}
	std::vector<std::vector<BackBoneSite>> helices(4);
	helices[0]=BackBoneBuilder::buildhelixat(15,dhelixstarts[1],direction,false);
	helices[1]=BackBoneBuilder::buildhelixat(15,dhelixstarts[2],direction,true);
//	helices[2]=BackBoneBuilder::buildhelixat(17,uhelixstarts[2],direction,true);
//	helices[3]=BackBoneBuilder::buildhelixat(17,uhelixstarts[1],direction,false);
/*
	std::vector<BackBoneSite> cterm=BackBoneBuilder::buildforwardbackbone(2,
				helices[1].back(),std::vector<std::pair<int,int>>(),
				std::vector<std::pair<int,int>>(),std::set<int>());
	std::vector<std::shared_ptr<std::vector<BackBoneSite>>> ls0s1,ls1s2,ls2s3,ls3h0,lh0h1;*/
	std::vector<std::shared_ptr<std::vector<BackBoneSite>>> ls2h0,lh0s0,ls0s1,ls1h1,lh1s3;
	int ntry=0;
	do { ls2h0=BackBoneBuilder::buildlinkers(5,strands[2].back(),helices[0][0],
		std::vector<std::pair<int,int>> (),
		std::vector<std::pair<int,int>>(),std::set<int>());
		++ntry;
		std::cout <<"Try linking s2 h0 " << ntry<<std::endl;
	}	while(ls2h0.empty() && ntry<10);
	if(ls2h0.empty()){
		std::cout<<"Failed to link s2 and h0"<<std::endl;
		exit(0);
	}
	ntry=0;
	do { lh0s0=BackBoneBuilder::buildlinkers(5,helices[0].back(),strands[0][0],
		std::vector<std::pair<int,int>> (),
		std::vector<std::pair<int,int>>(),std::set<int>());
		++ntry;
		std::cout <<"Try linking h0 s0 " << ntry<<std::endl;
	}	while(lh0s0.empty() && ntry<10);
	if(lh0s0.empty()){
		std::cout<<"Failed to link h0 and s0"<<std::endl;
		exit(0);
	}
	ntry=0;
	do { ls0s1=BackBoneBuilder::buildlinkers(5,strands[0].back(),strands[1][0],
		std::vector<std::pair<int,int>> (),
		std::vector<std::pair<int,int>>(),std::set<int>());
		++ntry;
		std::cout <<"Try linking s0 s1 " << ntry<<std::endl;
	}	while(ls0s1.empty() && ntry<10);
	if(ls0s1.empty()){
		std::cout<<"Failed to link s0 and s1"<<std::endl;
		exit(0);
	}
	ntry=0;
	do { ls1h1=BackBoneBuilder::buildlinkers(5,strands[1].back(),helices[1][0],
		std::vector<std::pair<int,int>> (),
		std::vector<std::pair<int,int>>(),std::set<int>());
		++ntry;
		std::cout <<"Try linking s1 h1 " << ntry<<std::endl;
	}	while(ls1h1.empty() && ntry<10);
	if(ls1h1.empty()){
		std::cout<<"Failed to link s1 and h1"<<std::endl;
		exit(0);
	}
	ntry=0;
	do { lh1s3=BackBoneBuilder::buildlinkers(5,helices[1].back(),strands[3][0],
		std::vector<std::pair<int,int>> (),
		std::vector<std::pair<int,int>>(),std::set<int>());
		++ntry;
		std::cout <<"Try linking h1 s3 " << ntry<<std::endl;
	}	while(lh1s3.empty() && ntry<10);
	ntry=0;
	if(lh1s3.empty()){
		std::cout<<"Failed to link h1 and s3"<<std::endl;
		exit(0);
	}

	std::vector<BackBoneSite> chain;
	for(auto &s:nterm) chain.push_back(s);
	for(auto &s:strands[2]) chain.push_back(s);
	for(auto &s:*(ls2h0[0])) chain.push_back(s);
	for(auto &s:helices[0]) chain.push_back(s);
	for(auto &s:*(lh0s0[0])) chain.push_back(s);
	for(auto &s:strands[0]) chain.push_back(s);
	for(auto &s:*(ls0s1[0])) chain.push_back(s);
	for(auto &s:strands[1]) chain.push_back(s);
	for(auto &s:*(ls1h1[0])) chain.push_back(s);
	for(auto &s:helices[1]) chain.push_back(s);
	for(auto &s:*(lh1s3[0])) chain.push_back(s);
	for(auto &s:strands[3]) chain.push_back(s);
	for(auto &s:cterm) chain.push_back(s);
	std::ofstream ofs;
	std::string filename="builtchain.pdb";
	ofs.open(filename.c_str());
	writeSitesToPDB(ofs,chain);
	ofs.close();
}
