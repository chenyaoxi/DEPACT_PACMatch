/*
 * nnregressionmodel.cpp
 *
 *  Created on: 2017年8月9日
 *      Author: hyliu
 */

#include "dstl/nnregressionmodel.h"

using namespace NSPdstl;
std::vector<double> LinearTransform::operator ()(const std::vector<double> & x) const {
	std::vector<double> result;
	assert(w_[0].size()==x.size());
	for(int row=0; row<b_.size(); ++row){
		double out=b_[row];
		const std::vector<double> & wr=w_[row];
		for(int col=0;col<x.size();++col) {
			out += wr[col]*x[col];
		}
		result.push_back(out);
	}
	return result;
}
double L3Estimator::operator ()(const std::vector<double> &input) const {
	std::vector<double> input2=conn12_(input);
	for(auto &val:input2) val=sigmoid(val);
	return conn23_(input2)[0];
}
void L3Estimator::setup(std::istream &is){
	conn12_=readlineartransform(is);
	conn23_=readlineartransform(is);
}
double L3Estimator::operator()(const std::vector<double> &input,
		std::vector<double> &dydx) const{
	std::vector<double> input2=conn12_(input);
	std::vector<double>  dsx2dx2;
	for (auto &val:input2){
		double dsdv;
		val=sigmoid(val,&dsdv);
		dsx2dx2.push_back(dsdv);
	}
	const std::vector<std::vector<double>> & dydsx2=conn23_.w();
	dydx.assign(input.size(),0.0);
	for(int j=0;j<input2.size();++j){
		double d=dydsx2[0][j]*dsx2dx2[j];
		const std::vector<double> & dx2jdx1=conn12_.w()[j];
		for(int i=0;i<input.size();++i){
			dydx[i] += d*dx2jdx1[i];
		}
	}
	return conn23_(input2)[0];
}

std::vector<double> L3SoftMaxClassifier::softmax(const std::vector<double> & logits,std::vector<std::vector<double>> *dpdlogits){
	std::vector<double> res;
	double ptot=0.0;
	for(auto lp:logits){
		double p=exp(lp);
		ptot +=p;
		res.push_back(exp(lp));
	}
	for(auto &p:res) {
		p/=ptot;
	}
	if(dpdlogits){
		dpdlogits->clear();
		int nclass=logits.size();
		for(int i=0;i<nclass;++i){
			dpdlogits->push_back(std::vector<double>(nclass));
			std::vector<double> &dpdlp=dpdlogits->back();
			for(int j=0;j<nclass;++j) {
				dpdlp[j]=-res[i]*res[j];
				if(i==j) dpdlp[j] +=res[i];
			}
		}
	}
	return res;
}
void L3SoftMaxClassifier::setup(std::istream &is){
	conn12_=readlineartransform(is);
	conn23_=readlineartransform(is);
}
std::vector<double> L3SoftMaxClassifier::operator()(const std::vector<double> &input,
		std::vector<std::vector<double>> &dpdx) const {
	std::vector<double> input2=conn12_(input);
	std::vector<double>  dsx2dx2;
	for (auto &val:input2){
		double dsdv;
		val=sigmoid(val,&dsdv);
		dsx2dx2.push_back(dsdv);
	}
	std::vector<double> logits=conn23_(input2);
	const std::vector<std::vector<double>> & dlpdsx2=conn23_.w();
	std::vector<std::vector<double>> dpdlps;
	std::vector<double> pr=softmax(logits,&dpdlps);

	int nclasses=pr.size();
	dpdx.assign(nclasses,std::vector<double>(input.size(),0.0));
	for(int n=0;n<nclasses;++n){
		std::vector<double> & dpndx=dpdx[n];
		std::vector<double> & dpndlp=dpdlps[n];
		for(int l=0;l<nclasses;++l){
			for(int j=0;j<input2.size();++j){
				double d=dpndlp[l]*dlpdsx2[l][j]*dsx2dx2[j];
				const std::vector<double> & dx2jdx1=conn12_.w()[j];
				for(int i=0;i<input.size();++i){
					dpndx[i] += d*dx2jdx1[i];
				}
			}
		}
	}
	return pr;
}

static int readintline(std::istream &is){
	char buffer[120];
	is.getline(buffer,120);
	return std::stoi(std::string(buffer));
}
static double readdoubleline(std::istream &is){
	char buffer[120];
	is.getline(buffer,120);
	return std::stod(std::string(buffer));
}
void NSPdstl::readvector(std::istream &is,std::vector<double> *vec){
	int size=readintline(is);
	vec->resize(size);
	for(double &val:*vec){
		val=readdoubleline(is);
	}
}
double CombinedEstimator::operator()(const std::vector<double> &input,std::vector<double> &dydx) const {
	double res=0;
	dydx.assign(input.size(),0.0);
	for(int i=0;i<estimators_.size();++i){
		std::vector<double> dydx_i;
		res +=estimators_[i](input,dydx_i);
		for(int l=0;l<dydx.size();++l){
			dydx[l]+=dydx_i[l];
		}
	}
	for(auto &v:dydx) v/=(double)estimators_.size();
	return res/(double) estimators_.size();
}

void CombinedEstimator::setup(std::istream &is){
	int nmodels=readintline(is);
	for(int i=0;i<nmodels;++i){
		estimators_.push_back(L3Estimator());
		estimators_.back().setup(is);
	}
}

LinearTransform NSPdstl::readlineartransform(std::istream &is){
	std::vector<std::vector<double>> w;
	std::vector<double> b;
	int nrows=readintline(is);
	for(int i=0;i<nrows;++i){
		w.push_back(std::vector<double>());
		readvector(is,&(w.back()));
	}
	readvector(is,&b);
	LinearTransform res;
	res.setup(w,b);
	return res;
}
L3Estimator NSPdstl::readl3estimator(std::istream &is){
	L3Estimator result;
	result.setup(is);
	return result;
}
L3Estimator NSPdstl::readl3estimator(const std::string &filename){
	std::ifstream ifs;
	ifs.open(filename.c_str());
	return readl3estimator(ifs);
	ifs.close();
}
L3SoftMaxClassifier NSPdstl::readl3softmaxclassifier(std::istream &is){
	L3SoftMaxClassifier result;
	result.setup(is);
	return result;
}

