/*
 * nnregressionmodel.h
 *
 *  Created on: 2017年8月9日
 *      Author: hyliu
 */

#ifndef DSTL_NNREGRESSIONMODEL_H_
#define DSTL_NNREGRESSIONMODEL_H_
#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

namespace NSPdstl{

inline double sigmoid(double x,double *dfdx=nullptr) {
	if(dfdx) *dfdx=0.0;
	if(x>1.e5) return 1.0;
	else if(x<-1.e5) return 0.0;
	double tmp=exp(x);
	double tmp1=tmp/(1.0+tmp);
	if(dfdx) {
		*dfdx=tmp1-tmp1*tmp1;
	}
	return tmp1;
};

inline double sigmoid(double x,double x0,double sigma,double *dfdx=nullptr){
	double res=sigmoid((x-x0)/sigma,dfdx);
	if(dfdx) *dfdx /=sigma;
	return res;
}
void readvector(std::istream &is,std::vector<double> *vec);
class LinearTransform {
public:
	std::vector<double> operator()(const std::vector<double> & x) const;
	void setup(const std::vector<std::vector<double>> & w, const std::vector<double> & b){
		assert(w.size()==b.size());
		w_=w;
		b_=b;
	}
	/**
	 *y=w_*x+b; dydx[i][j]=dy_i/dx_j=w_[i][j]
	 */
	const std::vector<std::vector<double>> & w() const {return w_;}
private:
	std::vector<std::vector<double>> w_;
	std::vector<double> b_;
};


class L3Estimator {
public:
	double operator()(const std::vector<double> &input) const;
	double operator()(const std::vector<double> &input,std::vector<double> &dydx) const;
	void setup(std::istream &is);

private:
	LinearTransform conn12_;
	LinearTransform conn23_;
};

class CombinedEstimator{
public:
	double operator()(const std::vector<double> &input) const {
		double res=0;
		for(auto & l3e:estimators_) res+=l3e(input);
		return res/(double) estimators_.size();
	}
	double operator()(const std::vector<double> &input,std::vector<double> &dydx) const;
	void setup(std::istream &is);
private:
	std::vector<L3Estimator> estimators_;
};

class L3SoftMaxClassifier {
public:
	std::vector<double> operator()(const std::vector<double> &input) const;
	std::vector<double> operator()(const std::vector<double> &input,
			std::vector<std::vector<double>> &dpdx) const ;
	void setup(std::istream &is);
	static std::vector<double> softmax(const std::vector<double> & logits,
			std::vector<std::vector<double>> *dpdlogits=nullptr);
private:
	LinearTransform conn12_;
	LinearTransform conn23_;
};

LinearTransform readlineartransform(std::istream &is);
L3Estimator readl3estimator(std::istream &is);
L3Estimator readl3estimator(const std::string &filename);
L3SoftMaxClassifier readl3softmaxclassifier(std::istream &is);
}




#endif /* DSTL_NNREGRESSIONMODEL_H_ */
