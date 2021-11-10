/*
 * pca.cpp
 *
 *  Created on: 2017年6月14日
 *      Author: hyliu
 */




#include "dstl/pca.h"
#include <cmath>

using namespace NSPdstl;
PCA::PCA(int ndim, bool angular) :
		ndim_(ndim), angular_(angular) {
	center_.resize(ndim, 0);
	covarmatrix_.resize(ndim, ndim, 0);
	components_.resize(ndim, ndim, 0);
	eigenvalues_.resize(ndim, 0);
}
std::vector<double> PCA::getPCAcrds(const std::vector<double> &point,int ndim) const {
	assert(ndim<=ndim_);
	std::vector<double> res;
	if(ndim==0) ndim=ndim_;
	for(int i=0; i<ndim; ++i) {
		const std::vector<double> & cpnti=getcomponent(i);
		double p{0};
		for(int j=0; j<ndim_;++j) {
			double delta= diff(point[j],center_[j]);
			p +=delta*cpnti[j];
		}
		res.push_back(p);
	}
	return res;
}
void PCA::calccenter(const std::vector<std::vector<double>> &data) {
	double datasize = 1.0 / (double) data.size();
	center_.resize(ndim_, 0.0);
	if (!angular_) {
		for (auto & d : data) {
			for (int j = 0; j < ndim_; j++) {
				center_[j] += d.at(j);
			}
		}
		for (auto &c : center_)
			c = c * datasize;
	} else {
		std::vector<double> xav(ndim_, 0.0);
		std::vector<double> yav(ndim_, 0.0);
		double deg = 3.1415926535897932384626 / 180.0;
		for (auto &d : data) {
			for (int j = 0; j < ndim_; j++) {
				double ang = d.at(j) * deg;
				xav[j] += cos(ang);
				yav[j] += sin(ang);
			}
		}
		for (int j = 0; j < ndim_; j++) {
			double x = xav[j] * datasize;
			double y = yav[j] * datasize;
			if (x > -1.e-10 && x < 1.e-10) {
				if (y >= 0)
					center_[j] = 90.0;
				else
					center_[j] = -90.0;
			} else {
				double ang = atan(y / x) / deg;
				if (x < 0)
					ang = ang + 180.0;
				if (ang > 180.0)
					ang -= 360.0;
				center_[j] = ang;
			}
		}
	}
}
void PCA::calccovarmatrix(const std::vector<double> &center,const std::vector<std::vector<double>> &data){
	center_=center;
	MyMatrix sum2dev;
	sum2dev.resize(ndim_,ndim_,0.0);
	for(auto &d:data) {
		for(int i=0; i<ndim_; ++i) {
			double devi=diff(d.at(i),center.at(i));
			sum2dev(i,i) += devi*devi;
			for(int j=0;j<i; ++j) {
				sum2dev(i,j) += devi*diff(d.at(j),center.at(j));
			}
		}
	}
	double ndata_i=1.0/(double) data.size();
	for(int i=0;i<ndim_; ++i) {
		for(int j=0; j<=i; ++j) {
			covarmatrix_(i,j)=sum2dev(i,j)*ndata_i;
			if(j !=i) covarmatrix_(j,i)=covarmatrix_(i,j);
		}
	}
}
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
MyMatrix NSPdstl::diagonizematrix(const MyMatrix &mat,std::vector<double> *eigenvalues){
		int ndim=mat.M;
	   gsl_matrix * matrix = gsl_matrix_alloc (ndim, ndim);
	   for(int i=0;i<ndim;++i)
		   for(int j=0; j<ndim;++j)
			   gsl_matrix_set(matrix,i,j,mat(i,j));
	   gsl_vector * eval = gsl_vector_alloc (ndim);
	   gsl_matrix * evec = gsl_matrix_alloc(ndim,ndim);
	   gsl_eigen_symmv_workspace * w=gsl_eigen_symmv_alloc (ndim);
	   gsl_eigen_symmv (matrix, eval,evec, w);
	   gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);
	   MyMatrix res;
	   res.resize(ndim,ndim,0.0);
	   eigenvalues->resize(ndim,0);
	   for (int i=0;i<ndim;i++) {
	 	   gsl_vector_view evec_i = gsl_matrix_column (evec, i);
	 	   std::vector<double> v;
	 	   for (int j=0;j<ndim;j++)
	 		     v.push_back(gsl_vector_get(&(evec_i.vector),j));
	   	   res.setcol(i,v);
	   	   (*eigenvalues)[i]=gsl_vector_get(eval,i);
	 	}
	   gsl_matrix_free(matrix);
	   gsl_matrix_free(evec);
	   gsl_vector_free(eval);
	   gsl_eigen_symmv_free(w);
	   return res;
}

