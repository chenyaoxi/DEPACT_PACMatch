/*
 * linkerproducer.h
 *
 *  Created on: 2016年11月24日
 *      Author: hyliu
 */

#ifndef BACKBONE_LINKERVENDOR_H_
#define BACKBONE_LINKERVENDOR_H_
#include <backbone/backbonelinker.h>
#include <backbone/backboneloop.h>
#include <map>
namespace NSPproteinrep {
class LinkerVendor {
public:
	LinkerVendor(){;}
	bool makeLinkers(unsigned int linkerkey, std::vector<std::vector<BackBoneSite>> *linker, unsigned int ncandidates,
			unsigned int ncopies=1,
			bool regenerate=false);
	 unsigned int bookloops(std::vector<std::vector<BackBoneSite>> *contextsegments,
			unsigned int headsegment, unsigned int tailsegment);
	bool makeLoops(unsigned int linkerkey, std::vector<std::vector<BackBoneSite>> *loop, unsigned int ncandidates,
				unsigned int ncopies=1,bool regenerate=false);
	unsigned int SegmentKey(const std::vector<BackBoneSite> *sg);
	unsigned int NewSegmentKey(const std::vector<BackBoneSite> *sg);
	unsigned int LinkerKey(unsigned int sg1key, unsigned int sg2key, unsigned int length);
private:
	static unsigned int MAXSEGMENTS;
	static unsigned int MAXCOPIES;
	std::map<const std::vector<BackBoneSite> *,unsigned int> segmentkeys_;
	std::vector<std::vector<BackBoneSite>> segmentheads_;
	std::vector<std::vector<BackBoneSite>> segmenttails_;
	std::map<unsigned int,BackBoneLinker> linkers_;
	std::map<unsigned int,BackBoneLoop> loops_;
	std::map<unsigned int, std::vector<std::vector<BackBoneSite>>> linkerstorage_;
	std::map<unsigned int, std::vector<std::vector<BackBoneSite>>> loopstorage_;
	void eraseLinkers(unsigned int segkey);
	static unsigned int genlinkerkey(unsigned int sg1key, unsigned int sg2key, unsigned int length) {
		return sg1key+sg2key*MAXSEGMENTS+length*MAXSEGMENTS*MAXSEGMENTS;
	}
	static unsigned int sgkey1(unsigned int linkerkey) {
		return linkerkey%MAXSEGMENTS;
	}
	static unsigned int sgkey2(unsigned int linkerkey) {
		return (linkerkey/MAXSEGMENTS)%MAXSEGMENTS;
	}
	static unsigned int linkerlength(unsigned int linkerkey){
		return linkerkey/(MAXSEGMENTS*MAXSEGMENTS);
	}
	void insertLinker(unsigned int linkerkey);
	void insertLoop(unsigned int linkerkey);
};
}






#endif /* BACKBONE_LINKERVENDOR_H_ */
