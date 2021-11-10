/*
 * segments.cpp
 *
 *  Created on: 2016年11月29日
 *      Author: hyliu
 */

#include <backbone/segments.h>
#include "dataio/splitstring.h"

using namespace NSPproteinrep;

void NSPproteinrep::readsegments(std::istream &ifs, std::vector<std::vector<BackBoneSite>> & segments){
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

