/*
 * segments.h
 *
 *  Created on: 2016年11月29日
 *      Author: hyliu
 */

#ifndef BACKBONE_SEGMENTS_H_
#define BACKBONE_SEGMENTS_H_
#include <backbone/backbonesite.h>

namespace NSPproteinrep {
void readsegments(std::istream &ifs, std::vector<std::vector<BackBoneSite>> & segments);
}



#endif /* BACKBONE_SEGMENTS_H_ */
