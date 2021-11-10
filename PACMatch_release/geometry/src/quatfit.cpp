/*
 * quatfit.cpp
 *
 *  Created on: 2016年4月12日
 *      Author: hyliu
 */

#include "geometry/quatfit.h"
#include "geometry/rotation.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>


using namespace NSPgeometry;
void QuatFit::transform(std::vector<XYZ> & crd){
	for(auto & p:crd) p=p+shift_;
	Rotation r;
//	QuaternionCrd q=getquaternioncrd();
//	std::cout <<q.angle() <<std::endl;
//	std::cout <<shift_.toString() <<std::endl;
	r.init(getquaternioncrd(),ref_center_);

//	RigidTransform rt(q,ref_center_,r.applytoCopy(shift_));
	for(auto & p:crd) {
		r.apply(&p);
	}
}
double QuatFit::fitting(const std::vector<double> & ref_coord,
						std::vector<double> & coord,std::vector<double> w){
	std::vector<XYZ> refcrd,crd;
	XYZ::vectortoxyzs(ref_coord,refcrd);
	XYZ::vectortoxyzs(coord,crd);
	double rmsd2=fitting(refcrd,crd,w);
	coord.clear();
	XYZ::xyzstovector(crd,coord);
	return rmsd2;
}

double QuatFit::fitting(const std::vector<XYZ> & ref_coord,
						std::vector<XYZ> & coord,std::vector<double> w){
	double rmsd2=setup(ref_coord, coord, w);
	transform(coord);
	return rmsd2;
}

QuaternionCrd QuatFit::getquaternioncrd(){
	return QuaternionCrd({m_evec[0][0],m_evec[1][0],m_evec[2][0],m_evec[3][0]});
}

double QuatFit::setup(const std::vector<double> &ref_coord,const std::vector<double> & coord,
		std::vector<double> w){
	std::vector<XYZ> refcrd, crd;
	XYZ::vectortoxyzs(ref_coord,refcrd);
	XYZ::vectortoxyzs(coord,crd);
	return setup(refcrd,crd,w);
}

double QuatFit::setup( std::vector<XYZ> ref_coord,
						std::vector<XYZ> coord,std::vector<double> w) {
	// The rotation that fit coord onto refcoord
	    gsl_matrix * B = gsl_matrix_alloc (4, 4);
	    gsl_matrix_set_zero (B);
	    double a[3],b[3];
	    for (unsigned int i=0;i<4;i++){
	    	for(unsigned int j=0;j<=i;j++){
	    		m_matrix[i][j]=0.0;
	    	}
	    }
	    int np=ref_coord.size();
	    if(w.empty()) {
	    	for(int i=0; i<np; ++i) w.push_back(1.0);
	    }
	    ref_center_=center(ref_coord,w);
//	    std::cout <<ref_center_.toString() <<std::endl;

	    shift_=-1.0*center(coord,w);

//	    std::cout <<shift_.toString() <<std::endl;
	    for(int i=0;i<np; ++i) {
	    	ref_coord[i] = ref_coord[i]-ref_center_;
	    	coord[i]=coord[i]+shift_;
	    }
	    shift_=shift_+ref_center_;

	    m_wtot=0.0;
	    int n=w.size();
	    for (int i=0; i<n;i++ ){
	    	m_wtot += w[i];
	    	for (int j=0; j<3;j++) {
	    		a[j]= ref_coord[i][j]+coord[i][j];
	    		b[j]= ref_coord[i][j]-coord[i][j];
	    	}
	    	m_matrix[0][0]+=w[i]*(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
	    	m_matrix[1][0]+=w[i]*(a[2]*b[1]-a[1]*b[2]);
	    	m_matrix[2][0]+=w[i]*(-a[2]*b[0]+a[0]*b[2]);
	    	m_matrix[3][0]+=w[i]*(a[1]*b[0]-a[0]*b[1]);

	    	m_matrix[1][1]+=w[i]*(b[0]*b[0]+a[2]*a[2]+a[1]*a[1]);
	    	m_matrix[2][1]+=w[i]*(b[1]*b[0]-a[0]*a[1]);
	    	m_matrix[3][1]+=w[i]*(b[2]*b[0]-a[0]*a[2]);

	    	m_matrix[2][2]+=w[i]*(b[1]*b[1]+a[2]*a[2]+a[0]*a[0]);
	    	m_matrix[3][2]+=w[i]*(b[2]*b[1]-a[1]*a[2]);

	    	m_matrix[3][3]+=w[i]*(b[2]*b[2]+a[1]*a[1]+a[0]*a[0]);
	    }
	    for (unsigned int i=0;i<4;i++){
	        double e=m_matrix[i][i]/m_wtot;
	        gsl_matrix_set(B,i,i,e);
	    	for(unsigned int j=0;j<i;j++){
	    		e=m_matrix[i][j]/m_wtot;
	    		gsl_matrix_set(B,i,j,e);
	    		gsl_matrix_set(B,j,i,e);
	    		m_matrix[j][i]=e;
	    	}
		}
	    gsl_vector *eval = gsl_vector_alloc (4);
	    gsl_matrix *evec = gsl_matrix_alloc (4,4);

	    gsl_eigen_symmv_workspace * ws = gsl_eigen_symmv_alloc (4);

	    gsl_eigen_symmv (B, eval, evec, ws);
	    gsl_eigen_symmv_free(ws);

	    gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);
	    for(unsigned int i=0;i<4;i++){
	    	m_eval[i]=gsl_vector_get(eval,i);
	    	for(unsigned int j=0;j<4;j++)
	    		m_evec[i][j]=gsl_matrix_get(evec,i,j);
	    }
	    gsl_matrix_free (B);
	    gsl_matrix_free (evec);
	    gsl_vector_free (eval);
	    return m_eval[0];
}


