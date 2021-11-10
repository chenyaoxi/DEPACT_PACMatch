/*
 * clusterconformers.cpp
 *
 *  Created on: 2017年6月28日
 *      Author: hyliu
 */
#include "dataio/pdbfolder.h"
#include "dataio/inputlines.h"
#include "proteinrep/aaconformer.h"
#include "kmeans/runkmeans.h"
#include <string>
using namespace NSPdataio;
using namespace NSPproteinrep;

int main(int argc, char **argv){
	TextLines txts;
	txts.init(std::string(argv[1]));
	std::multimap<std::string,char> pdbchainids;
	for(auto & word:txts.lines()){
		pdbchainids.insert(
				std::make_pair(word.substr(0,4),word[4]));
	}
	std::string pdbroot="";
	if(argc>2)pdbroot=std::string(argv[2]);
	PdbFolder &pdbfolder= PdbFolder::getinstance(pdbroot);
	std::string selectresidue="GLU";
	std::vector<AAConformer> selectedconformers;
	for(auto & pdbchain:pdbchainids){
		std::string pdbid=pdbchain.first;
		char chainid=pdbchain.second;
		std::cout <<"Reading " <<pdbid << " chain " <<chainid<<std::endl;
		std::string localfilename=pdbfolder.getfile(pdbid);
		if(localfilename.empty()) {
			std::cout <<"Get PDB file for " << pdbid <<"failed." <<std::endl;
			abort();
		}
		PDBModelID mid;
		mid.pdbid=pdbid;
		AAConformersInModel modelconformers;
		modelconformers.pdbmodelid=mid;
		modelconformers.readpdbfile(localfilename);
		if (chainid=='A') chainid=' ';
		const std::vector<AAConformer> & conformers=modelconformers.conformersinchain(chainid);
		for(auto & c:conformers) {
			if(c.residuename != selectresidue) continue;
			if(!c.sidechaincrdcomplete()) continue;
			selectedconformers.push_back(c);
			selectedconformers.back().removeoxtorother("O");
			selectedconformers.back().removeatomsnotlisted();
		}
	}
	NSPkmeans::RunKMeans rkmeans;
	int dim=3*selectedconformers[0].getglobalcrd().size();
	int maxpoints=selectedconformers.size();
	rkmeans.init_data(dim,maxpoints);
	for(auto &c:selectedconformers){
		std::map<std::string,NSPgeometry::XYZ> & localcrd=c.calclocalcrd();
		std::vector<double> point;
		if(localcrd.size() != dim/3) {
			std::cout <<"crd size "<<localcrd.size()<<std::endl;
		}
		for(auto &atm:localcrd) {
			point.push_back(atm.second.x_);
			point.push_back(atm.second.y_);
			point.push_back(atm.second.z_);
		}
		rkmeans.addpoint(point);
	}
	rkmeans.finishdata();
	for(int k=3; k<40; ++k) {
		rkmeans.setrun(k,NSPkmeans::RunKMeans::HYBRID);
		rkmeans.run();
		std::cout<<k<<'\t'<<rkmeans.getdistortion()<<std::endl;
/*	auto &clusters=rkmeans.clusters();
	for(auto &c:clusters){
		std::cout<<c.distortion;
		std::cout<<'\t'<<c.members.size()<<std::endl;
		for(int i=0; i<dim;++i)
		std::cout <<'\t'<< c.center[i];
		std::cout<<std::endl;
	}*/
	}
}
