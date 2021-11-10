/*
 * backbone.h
 *
 *  Created on: 2016年11月17日
 *      Author: hyliu
 */

#ifndef BACKBONE_BACKBONEDATA_H_
#define BACKBONE_BACKBONEDATA_H_
#include "proteinrep/idealgeometries.h"
namespace NSPproteinrep {

class BackBoneData {
public:
	static IdealGeometries & idealgeometries(const std::string & filename=" "){
		static IdealGeometries instance_;
		if(filename !=" " && filename != instance_.filename()) {
			instance_.readIdealValues(filename);
		}
		return instance_;
	}
	static std::map<int,std::string> atomnames;// initialize to {"N","CA","C","O","H"};
	enum {N=0,CA=1,C=2,O=3,H=4};
private:
};

}




#endif /* BACKBONE_BACKBONEDATA_H_ */
