/*
 * testlinkervendor.cpp
 *
 *  Created on: 2016年11月24日
 *      Author: hyliu
 */

#include <backbone/linkervendor.h>
#include "dataio/splitstring.h"
using namespace NSPproteinrep;
void readsegments(std::istream &ifs, std::vector<std::vector<BackBoneSite>> & segments){
	char buffer[120];
	ifs.getline(buffer,120);
	std::vector<int> segmentlengths=NSPdataio::integersInString(std::string(buffer));
	segments.resize(segmentlengths.size(),std::vector<BackBoneSite>());
	for(int i=0;i<segmentlengths.size();++i) {
		if(!readbackbonesites(ifs,segmentlengths[i],segments.at(i))){
			std::cout <<"Error reading segments."<<std::endl;
			exit(1);
		}
	}
}
int main(int argc, char **argv){
	if(argc < 2) {
		std::cout <<"Usage: " <<argv[0] <<" segments_filename "<<std::endl;
		exit(1);
	}
	IdealGeometries & idg=IdealGeometries::getGlobalInstance("idealgeometries.dat");
	std::vector<std::vector<BackBoneSite>> segments;
	std::ifstream ifs;
	ifs.open(argv[1]);
	readsegments(ifs,segments);
	ifs.close();
	LinkerVendor lv;
	for(auto it=segments.begin(); it != segments.end(); ++it) {
		lv.SegmentKey( &(*it));
	}
// test genlinker
	std::vector<std::vector<BackBoneSite>> loops;
	unsigned int ncandidates=500;
	unsigned int ncopies=10;
	unsigned int k1=lv.SegmentKey(&segments[0]);
	unsigned int k2=lv.SegmentKey(&segments[1]);
	unsigned int linkerkey=lv.LinkerKey(k1,k2,5);
	bool success=lv.makeLoops(linkerkey,
			&loops,ncandidates,ncopies,
			false);
	if(!success) {
		std::cout <<"make Linker failed" <<std::endl;
		exit(0);
	}
// write pdb files
	std::ofstream ofs;
	for(int i=0; i<ncopies; ++i){
		std::string filename="sg1-loop"+std::to_string(i)+"sg2.pdb";
		ofs.open(filename.c_str());
		std::vector<BackBoneSite> collection;
		collection.resize(segments[0].size()+loops[i].size()+segments[1].size()-2);
		std::copy(segments[0].begin(),segments[0].end()-1,collection.begin());
		std::copy(loops[i].begin(),loops[i].end(),collection.begin()+segments[0].size()-1);
		std::copy(segments[1].begin()+1,segments[1].end(),collection.begin()+segments[0].size()+loops[i].size()-1);
		writeSitesToPDB(ofs,collection);
		ofs.close();
	}
}


