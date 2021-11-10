/*
 * linkerbuildercontrol.h
 *
 *  Created on: 2016年11月29日
 *      Author: hyliu
 */

#ifndef BACKBONE_LINKERBUILDERCONTROLS_H_
#define BACKBONE_LINKERBUILDERCONTROLS_H_

#include "dataio/parameters.h"

namespace NSPproteinrep {

class LinkerBuilderControls {
public:
	static void initcontrols(const std::string &filename);
	static LinkerBuilderControls & getinstance() {
		static LinkerBuilderControls instance;
		return instance;
	}
	std::string segmentsfilename;
	std::string outlinkerfilename;
	std::string idealgeometriesfilename;
	unsigned int ncandidates;
	std::vector<std::vector<unsigned int>> lengthsandcopies;
private:
	LinkerBuilderControls(){;}
};

}

#endif /* BACKBONE_LINKERBUILDERCONTROLS_H_ */
