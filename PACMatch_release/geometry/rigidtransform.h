/*
 * rigidtranform.h
 *
 *  Created on: 2016年11月8日
 *      Author: hyliu
 */

#ifndef RIGIDTRANSFORM_H_
#define RIGIDTRANSFORM_H_

#include "geometry/rotation.h"

namespace NSPgeometry {
/**A rigidtransform comprises a rotation and a translation.
 *
 * X_transformed=rotation(X_original)+translation
 */
class RigidTransform{

public:

	RigidTransform(){;}
	/** constructor
	 *
	 * @param Q: the roation represented as a quternion
	 * @param rotcenter : a point on the rotation axis
	 * @param t the translation
	 */
	RigidTransform(const QuaternionCrd &Q, const XYZ &rotcenter=XYZ(), const XYZ &t=XYZ()):trans_(t){
		rot_.init(Q,rotcenter);
	}
	/** constructor
	 *
	 * @param r: the rotation
	 * @param t: the translation
	 */
	RigidTransform(const Rotation & r, const XYZ & t=XYZ()): rot_(r),trans_(t){;}
	/** apply the rigid transform to a point
	 *
	 * @param p: point to the XYZ point to be transformed
	 */
	void apply(XYZ *p) const;

	/**Apply the rigid transform to a copy of an input point
	 *
	 * @param p: input XYZ point, will not be changed.
	 * @return the transformed XYZ point.
	 */
	XYZ applytoCopy(const XYZ & p) const;
	/**get the rotation part of the rigid transform
	 *
	 */
	Rotation &rotation() {return rot_;}
	/** get the const rotation part of the rigid transform
	 *
	 */
	const Rotation & rotation() const {return rot_;}
	/** get the translation part of the rigid transform
	 *
	 */
	XYZ & translation() {return trans_;}
	/** get the const translation oart of the rigid transform
	 *
	 */
	const XYZ & translation() const {return trans_;}
private:
	/**the rotation.It co
	 *
	 */
	Rotation rot_;

	/**The translatio
	 *
	 */
	XYZ trans_;
};
/**Combine a rigid transform with a rotation, get a new rigid transform
 * X_transformed=Rotation_r(RigidTranform rt(X_original)
 * @return the resulting rigid transform
 */
RigidTransform applyRotation(const Rotation & r, const RigidTransform &rt);
/**Generate a rigid transformation with random roation and translation
 * The rotation is around a random axis passing the origin, with uniform distributed angles between
 * 0 and maxrotate(in degrees)
 * The translation is of length between 0 to maxtranslate
 *
 */
RigidTransform randomrigidtransform(double maxrotate, double maxtranslate);

/**Get the reigidtransformation that would superimpose the coordinate set crdb onto the
 * coordinate set crda.
 * @param crda: the reference(fixed) coordinate
 * @param crdb: the coordinates to be superimposed (only the transform is calculated, the
 * coorinates are not changed.
 * @param alignedpositions: the aligned positions between crda and crdb. Only the atoms in the aligned positions will be
 * considered for superposition.
 * @param dev2: if not null, the squared rmsd will be delivered to dev2.
 *
 */
RigidTransform superpose(const std::vector<XYZ> & crda, const std::vector<XYZ> &crdb,
		const std::vector<std::pair<int,int>> & alignedpostions,double *dev2=nullptr);
}



#endif /* RIGIDTRANSFORM_H_ */
