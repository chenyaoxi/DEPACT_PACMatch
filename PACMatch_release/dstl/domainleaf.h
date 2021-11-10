/*
 * makedomaintree.h
 *
 *  Created on: 2016年2月26日
 *      Author: hyliu
 */

#ifndef DOMAINLEAF_H_
#define DOMAINLEAF_H_
#include "dstl/domaintree.h"
namespace domaintree {
template<typename T = long>
struct LeafStruct {
	std::vector<T> points;
};

template<typename T = long>
class LeafOperation {
public:
	void operator()(Domain & d, std::shared_ptr<LeafStruct<T>> & leaf) {
		if (!leaf) {
			leaf = std::shared_ptr < LeafStruct < T >> (new LeafStruct<T>());
		}
		leaf->points.push_back(current_point);
//		d.print();
	}
	LeafOperation(T p) :
			current_point(p) {
		;
	}
private:
	T current_point;
};

struct UsualCrd {
	static double diff(double d1, double d2) {
		return d1 - d2;
	}
	static std::vector<double> diff(const std::vector<double> & d1, const std::vector<double> & d2) {
		std::vector<double> delta;
		for(int d=0; d<d1.size();++d) delta.push_back(d1[d]-d2[d]);
		return delta;
	}


	static std::vector<double> shift(std::vector<double> a,const std::vector<double> & base){
		return a;
	}

	static double mean(double sum, double count){
		sum=sum/count;
		return sum;
	}
	static std::vector<double> mean(std::vector<double> sum, double count){
		for(int d=0; d<sum.size(); ++d) {
			sum[d]=sum[d]/count;
		}
		return sum;
	}
};

struct AngleCrd {
	static double diff(double d1, double d2) {
		double tmp = d1 - d2;
		while (tmp < -180)
			tmp += 360;
		while (tmp > 180)
			tmp -= 360;
		return tmp;
	}
	static std::vector<double> diff(const std::vector<double> & d1, const std::vector<double> & d2) {
		std::vector<double> delta;
		for(int d=0; d<d1.size(); ++d) delta.push_back(diff(d1[d],d2[d]));
		return delta;
	}

	static double shift(double a,double base){
		if(a-base >180.0) a -=360.0;
		if(a-base <-180.0) a +=360.0;
		return a;
	}


	static std::vector<double> shift(std::vector<double> a,const std::vector<double> & base){
		for(int d=0; d<a.size();++d){
			if(a[d]-base[d] >180.0) a[d] -=360.0;
			if(a[d]-base[d] <-180.0) a[d] +=360.0;
		}
		return a;
	}

	static double mean(double sum, double count){
		sum=sum/count;
		if(sum >180.0) sum -=360.0;
		if(sum <-180.0) sum +=360.0;
		return sum;
	}
	static std::vector<double> mean(std::vector<double> sum, double count){
		for(int d=0; d<sum.size(); ++d) {
			sum[d]=sum[d]/count;
			if(sum[d] >180.0) sum[d] -=360.0;
			if(sum[d]<-180.0) sum[d] +=360.0;
		}
		return sum;
	}
	static std::vector<double> average(std::vector<std::vector<double>> &data,
			std::vector<double> &var){
		if(data.empty()) return std::vector<double>();
		int ndim=data[0].size();
		std::vector<double> sum(ndim,0.0);
		double count;
		for(long pidx=0; pidx <data.size();++pidx) {
				std::vector<double> tmp=shift(data[pidx],data[0]);
				for (int d=0; d<ndim; ++d)
					sum[d] += tmp[d];
				count +=1.0;
		}
		std::vector<double> av=mean(sum,count);
		var.resize(ndim,0.0);
		for (long i=0; i<data.size(); ++i) {
			std::vector<double> delta=diff(av,data[i]);
			for(int d=0; d<ndim; ++d) var[d] += delta[d]*delta[d];
		}
		for(int d=0; d<ndim; ++d) var[d] = var[d]/count;
		return av;
	}

};
template<typename T = long, typename DATA = std::vector<Domain::point_type>,
		typename CRDTYPE = UsualCrd,typename NN=NNearest<T>>
class D2Leaf {

public:
	typedef CRDTYPE crd_type;

//	D2Leaf(DATA *d) :
//			data_(d), nearest_neighbor_(-1){
//		;
//	}

	D2Leaf(DATA *d, int n=1, double cut2=100000000000.0):
		data_(d), nearest_neighbor_(-1),nnearest_(n,cut2){;}

	void setpoint(const Domain::point_type & p) {
		current_point_ = p;
	}
	T getneighbour() {
		return nearest_neighbor_;
	}
	NN & nnearest() {
		return nnearest_;
	}
	const NN &nnearest() const {
		return nnearest_;
	}
	void operator()(std::shared_ptr<LeafStruct<T>> & leaf, double & disbound) {
		if (disbound < 0)
			disbound = nnearest_.maxd2();
		for (auto i : leaf->points) {
			double dis2 = 0.0;
			std::vector<double> tmp=crd_type::diff((*data_)[i], current_point_);
			for (int d = 0; d < current_point_.size(); ++d) {
//				double tmp = crd_type::diff((*data_)[i][d], current_point_[d]);
				dis2 += tmp[d] * tmp[d];
			}
			if(nnearest_.insert(i,dis2,& disbound)){
				nearest_neighbor_=nnearest_.nearestone();
//				std::cout <<"disbound: " << disbound <<"\n";
			}
		}
	}
	void operator()(std::shared_ptr<LeafStruct<T>> & leaf, const Domain::point_type &point,
			double & disbound) {
		setpoint(point);
		(*this)(leaf,disbound);
	}
private:
	Domain::point_type current_point_;
	T nearest_neighbor_;
	NN nnearest_;
	DATA* data_;
};
}
#endif /* DOMAINLEAF_H_ */
