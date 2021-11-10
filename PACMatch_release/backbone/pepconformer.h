/*
 * pepconformer.h
 *
 *  Created on: 2017年6月29日
 *      Author: hyliu
 */

#ifndef BACKBONE_PEPCONFORMER_H_
#define BACKBONE_PEPCONFORMER_H_
#include "geometry/localframe.h"
#include "backbone/backbonesite.h"

namespace NSPproteinrep {
class PepConformer {
public:
	int & length(){return length_;}
	const int & length() const {return length_;}
	std::vector<NSPgeometry::XYZ> &globalcrd(){return globalcrd_;}
	std::vector<NSPgeometry::XYZ> &localcrd(){return localcrd_;}
	const std::vector<NSPgeometry::XYZ> &globalcrd() const {return globalcrd_;}
	const std::vector<NSPgeometry::XYZ> &localcrd() const {return localcrd_;}
	NSPgeometry::LocalFrame & localframe(){return localframe_;}
	const NSPgeometry::LocalFrame &localframe() const {return localframe_;}

	NSPgeometry::LocalFrame &calclocalframe();
	static NSPgeometry::LocalFrame calclocalframe(
			const NSPgeometry::XYZ & rca,const NSPgeometry::XYZ & rn, const NSPgeometry::XYZ & rc);
	std::vector<NSPgeometry::XYZ> & calclocalcrd();
	void tobackbonesites(std::vector<BackBoneSite> &sites,
			const std::vector<NSPgeometry::XYZ> & crd);
private:
	int length_;
	std::vector<NSPgeometry::XYZ> globalcrd_;
	std::vector<NSPgeometry::XYZ> localcrd_;
	NSPgeometry::LocalFrame localframe_;
};

PepConformer make_conformer(int length,const std::vector<BackBoneSite> & sites, int posi);
}




#endif /* BACKBONE_PEPCONFORMER_H_ */
