/*
 * pca.h
 *
 *  Created on: 2017年6月14日
 *      Author: hyliu
 */

#ifndef DSTL_PCA_H_
#define DSTL_PCA_H_

#include <vector>
#include <algorithm>
#include <cassert>
namespace NSPdstl {
template<typename T = double>
struct MatrixAsColvec {
	int M { 0 }, N { 0 };
	std::vector<std::vector<T>> vectors;
	T & operator()(int i, int j) {
		return vectors.at(j).at(i);
	}
	const T & operator()(int i, int j) const {
		return vectors.at(j).at(i);
	}
	MatrixAsColvec() {
		;
	}
	MatrixAsColvec(int dim1, int dim2, const T & ele = T()) :
			M(dim1), N(dim2) {
		vectors.resize(dim2, std::vector<T>(dim1, ele));
		;
	}
	void resize(int dim1, int dim2, const T & ele = T()) {
		M=dim1;
		N=dim2;
		vectors.resize(dim2);
		for (auto &v : vectors)
			v.resize(dim1, ele);
	}
	void setrow(int i, const std::vector<T> & inrow) {
		assert(i < M);
		int j = 0;
		for (auto &v : vectors)
			v[i] = inrow[j++];
	}
	void setcol(int j, const std::vector<T> & incol) {
		assert(j < N);
		assert(incol.size() == M);
		std::copy(incol.begin(), incol.end(), vectors[j].begin());
	}
	std::vector<T> getrow(int i) const {
		assert(i < M);
		std::vector<T> row;
		for (auto &v : vectors)
			row.push_back(v[i]);
		return row;
	}
	std::vector<T> & getcol(int j) {
		assert(j < N);
		return vectors[j];
	}
	const std::vector<T> & getcol(int j) const {
		assert(j < N);
		return vectors[j];
	}
	MatrixAsColvec<T> operator *(const MatrixAsColvec<T> mat2) {
		assert(N == mat2.M);
		MatrixAsColvec<T> result(M, mat2.N, 0);
		for (int j = 0; j < mat2.N; ++j) {
			const std::vector<T> & m2colj = mat2.getcol(j);
			std::vector<T> & rescolj = result.getcol(j);
			for (int i = 0; i < M; ++i) {
				T sum { 0 };
				for (int k = 0; k < N; ++k) {
					sum += vectors[k][i] * m2colj[k];
				}
				rescolj[i] = sum;
			}
		}
		return result;
	}
	void transpose() {
		(*this) = this->gettranspose();
	}
	MatrixAsColvec<T> gettranspose() const {
		MatrixAsColvec<T> result(N, M);
		int j = 0;
		for (auto & v : vectors)
			result.setrow(j++, v);
		return result;
	}
};
typedef MatrixAsColvec<double> MyMatrix;
MyMatrix diagonizematrix(const MyMatrix &mat, std::vector<double> *eigenvalues);
class PCA {
public:
	PCA(int ndim, bool angular = false);

	void doPCA(const std::vector<std::vector<double>> &data) {
		calccovarmatrix(data);
		doPCA();
	}
	void doPCA(const std::vector<double> &center,
			const std::vector<std::vector<double>> &data) {
		calccovarmatrix(center, data);
		doPCA();
	}

	void doPCA() {
		components_ = diagonizematrix(covarmatrix_, &eigenvalues_);
	}
	void calccenter(const std::vector<std::vector<double>> &data);
	void calccovarmatrix(const std::vector<std::vector<double>> &data) {
		calccenter(data);
		calccovarmatrix(center_, data);
	}
	void calccovarmatrix(const std::vector<double> &center,
			const std::vector<std::vector<double>> &data);
	const std::vector<double> &getcenter() const {
		return center_;
	}
	double geteigenvalue(int i) const {
		return eigenvalues_[i];
	}
	const std::vector<double> & getcomponent(int i) const {
		return components_.getcol(i);
	}
	double getcovar(int i,int j) {
		return covarmatrix_(i,j);
	}
	std::vector<double> getPCAcrds(const std::vector<double> &point, int ndim =
			0) const;
private:
	int ndim_;
	bool angular_ { false };
	std::vector<double> center_;
	MyMatrix covarmatrix_;
	MyMatrix components_;
	std::vector<double> eigenvalues_;
	double diff(double x1, double x2) const {
		if (!angular_)
			return x1 - x2;
		else {
			double d = x1 - x2;
			while (d > 180.0)
				d -= 360.0;
			while (d < -180.0)
				d += 360.0;
			return d;
		}
	}
};

}




#endif /* DSTL_PCA_H_ */
