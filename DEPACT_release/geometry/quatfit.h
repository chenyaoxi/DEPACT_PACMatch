/*
 * quatfit.h
 *
 *  Created on: 2016年4月12日
 *      Author: hyliu
 */

#ifndef QUATFIT_H_
#define QUATFIT_H_
#include "quaternioncrd.h"
#include "rigidtransform.h"
namespace NSPgeometry {
/*!the class carries out minimum RMSD fit between two sets of Cartesian
 * Coordinates using the quaternion approach.
 *
 */
class QuatFit {
public:
	QuatFit():m_wtot(0.0){}
	/*!calculate and store the transformation to fit coord to ref_coord.
	 * Actual coordinate transformation is not executed in this function.
	 * Input 3-d coordinates are packed as one-dimensional vectors
	 * weights are of one-third of the length these 1-d vectors.
	 * returns rmsd squared.
	 */
	double setup(const std::vector<double> &ref_coord,const std::vector<double> & coord,
			 std::vector<double> weights=std::vector<double>());
	/*!calculate and store the transformation to fit coord to ref_coord.
	 * Actual coordinate transformation is not executed in this function.
	 * Input 3-d coordinates are vectors of XYZ objects
	 * returns rmsd squared.
	 *
	 */
	double setup(std::vector<XYZ> ref_coord,std::vector<XYZ> coord,
			const std::vector<double> weights=std::vector<double>());

	double evec(int i,int j) {return m_evec[i][j];}
	double matrix(int i, int j) {return m_matrix[i][j];}
	/*!
	 * get the rotation part
	 */
	QuaternionCrd getquaternioncrd();

	/*!
	 * calculate and store minimum RMSD transformations, apply them to coord.
	 * 3-d coordinates are packed into 1-d vectors.
	 * weights are one-third of the length of these 1-d vectors
	 * returns rmsd squared.
	 */
	double fitting(const std::vector<double> &ref_coord,std::vector<double> & coord,
			std::vector<double> weights=std::vector<double>());
	/*!
	 * calculate and store minimum RMSD transformations, apply them to coord.
	 * returns rmsd squared.
	 */
	double fitting(const std::vector<XYZ> &ref_coord,std::vector<XYZ> & coord,
			std::vector<double> weights=std::vector<double>());

	/*!
	 * apply the stored transformation to crd
	 */
	void transform(std::vector<XYZ> & crd);

	/*!
	 * get the transformation as a RigidTransform object.
	 */
	RigidTransform getRigidTransform(){
		QuaternionCrd Q=getquaternioncrd();
		Rotation r0(Q,XYZ());
		Rotation r(Q,ref_center_);
		RigidTransform rt(Q,ref_center_,r0.applytoCopy(shift_));
		return rt;
	}
private:
	/*!eigen value part of the minimum RMSD solution
	 * m_eval[0] is the minimum rmsd squared.
	 *
	 */
	double m_eval[4];

	/*!
	 * eigen vector part of the minimum RMSD solution
	 * m_evec[0] is the rotation quaternion
	 */
	double m_evec[4][4];
	/*!
	 * rotation matrix
	 */
	double m_matrix[4][4];
	/*!
	 * total weight
	 */
	double m_wtot;
	/*!
	 * (weighted) geometric center of the reference coordinates
	 */
	XYZ ref_center_;
	/*!
	 * translation part of the fit
	 */
	XYZ shift_;
};
}

#endif /* QUATFIT_H_ */
