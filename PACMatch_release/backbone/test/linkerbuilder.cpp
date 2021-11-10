/*
 * likerbuilder.cpp
 *
 *  Created on: 2016年11月29日
 *      Author: hyliu
 */
#include "backbone/segments.h"
#include "backbone/linkerbuildercontrols.h"
#include "backbone/linkervendor.h"
#include <iostream>
#include <fstream>

using namespace NSPproteinrep;
using namespace NSPdataio;
int main(int argc, char **argv) {
	LinkerBuilderControls::initcontrols(std::string(argv[1]));
	LinkerBuilderControls & par=LinkerBuilderControls::getinstance();

	IdealGeometries::getGlobalInstance(par.idealgeometriesfilename);

	std::ifstream ifs;
	ifs.open(par.segmentsfilename.c_str());
	std::vector<std::vector<BackBoneSite>> segments;
	readsegments(ifs,segments);
	ifs.close();

	LinkerVendor lv;
	std::ofstream ofs;
	ofs.open(par.outlinkerfilename.c_str());
//	std::ostream & ofs= std::cout;
	for (auto & lc: par.lengthsandcopies) {
		std::vector<std::vector<BackBoneSite>> linkers;
		unsigned int seg1=lc[0];
		unsigned int seg2=lc[1];
		unsigned int length=lc[2];
		unsigned int ncopies=lc[3];
		unsigned int linkerkey=lv.LinkerKey(lv.SegmentKey(&segments[seg1]), lv.SegmentKey(&segments[seg2]),length);
		bool success=lv.makeLinkers(linkerkey, &linkers,par.ncandidates,ncopies,true);
		if(!success) {
				std::cout <<"make Linker failed for segments " << seg1 <<" " <<seg2 <<std::endl;
		} else {
			ofs << seg1 <<" " <<" " <<seg2 << " " << length << " " <<ncopies <<std::endl;
			for(auto & l:linkers) {
				for( auto & s:l) {
					ofs <<s.toString();
				}
			}
		}
	}
	ofs.close();
}
