/*
 * dpclustering.h
 *
 *  Created on: 2019年10月10日
 *      Author: hyliu
 */

#ifndef DPCLUSTERING_H_
#define DPCLUSTERING_H_
#include <vector>
#include <cassert>
#include <iostream>
#include <math.h>
#include "dstl/sortindex.h"
namespace NSPdstl{
/* Template class for density-peak clustering
 * Reference:
 * Alex Rodriguez and Alessandro Laio, Clustering by fast search and find of
 * density peaks, Science 344, 1492(2014)
 * Each item to be clustered is of type T
 * To use the template class, the class must be able to call the following function
 * to calculate the distance between two points t1 and t2:
 * double distance(const T & t1, const T &t2)
 * If T correspond to type std::vector<double>, this distance function has been
 * internally defined
 *
 * The kernel function to estimate local density from the distance is
 * exp(-d^2/sigma^2), sigma is a parameter provided to the constructor
 *
 * See example codes in testdpclustering.cpp for how to use the class
 *
 */
template <typename T>
class DPClustering{
public:
	DPClustering (double sigma){sigma2_=sigma*sigma;}
	void setdistmatrix(const std::vector<T> & nodes){
		nnodes_=nodes.size();
		distmatrix_.clear();
		for(int i=0;i<nnodes_;++i){
			distmatrix_.push_back(std::vector<double>());
			for(int j=0;j<=i;++j){
				double d=distance(nodes.at(i),nodes.at(j));
				distmatrix_.back().push_back(d);
			}
		}
	}
	void setdistmatrix_hm(const std::vector<T> & nodes){
		nnodes_=nodes.size();
		distmatrix_.clear();
		for(int i=0;i<nnodes_;++i){
			distmatrix_.push_back(std::vector<double>());
			for(int j=0;j<=i;++j){
				double d=distance_hm(nodes.at(i),nodes.at(j));
				distmatrix_.back().push_back(d);
			}
		}
	}
	void  setdensityidx(){
		assert(!distmatrix_.empty());
		density_.assign(nnodes_,0.0);
		for(int i=1; i <nnodes_;++i ){
			for (auto j=0; j <i ;++j){
				double d=distmatrix_[i][j];
				double o=exp(-d*d/sigma2_);
				density_[i] -=o;
				density_[j] -=o;
			}
		}
		densityidx_=sortindex(density_);
		for(auto &d:density_) d=-d;
		dbetter_.resize(nnodes_);
		for(int i=1;i<nnodes_;++i){
			int idxi=densityidx_[i];
			double dmin=1000000000.0;
			int imin=-1;
			for(int j=0;j<i;++j){
				int idxj=densityidx_[j];
				double dij;
				if(idxi>idxj) dij=distmatrix_[idxi][idxj];
				else dij=distmatrix_[idxj][idxi];
				if(dij<dmin) {
					dmin=dij;
					imin=idxj;
				}
			}
			dbetter_[idxi]=std::make_pair(imin,dmin);
		}
	}
	const std::vector<int> & docluster(double dminbetweencenter){
		clusterid_.assign(nnodes_,-1);
		clusterid_[densityidx_[0]]=0;
		centerindices_.clear();
		centerindices_.push_back(densityidx_[0]);
		int maxid=0;
		for(int i=1;i<nnodes_;++i){
			int idxi=densityidx_[i];
			if(dbetter_[idxi].second > dminbetweencenter){
				clusterid_[idxi]=++maxid;
				centerindices_.push_back(idxi);
			} else {
				int idxj=dbetter_[idxi].first;
				assert(clusterid_[idxj]>=0);
				clusterid_[idxi]=clusterid_[idxj];
			}
		}
		nclusters_=maxid+1;
		return clusterid_;
	}
	const std::vector<int> & docluster(const std::vector<T> & nodes,double dminbetweencenter){
		setdistmatrix(nodes);
		setdensityidx();
		return docluster(dminbetweencenter);
	}
	const std::vector<int> & docluster_hm(const std::vector<T> & nodes,double dminbetweencenter){
		setdistmatrix_hm(nodes);
		setdensityidx();
		return docluster(dminbetweencenter);
	}
	const std::vector<int> & clusterid() const{return clusterid_;}
	int nclusters() const {return nclusters_;}
	const std::vector<int> &centers() const {return centerindices_;}
	void makeplot(std::ostream & os) const {
		for(int i=0;i<nnodes_;++i){
			os << density_.at(i)<<"   "<< dbetter_.at(i).first<<std::endl;
		}
	}
	double distance(const std::vector<double>& p1, const std::vector<double>& p2){
		double d=0.0;
		for(int i=0;i<p1.size();++i){
			double delta=p1[i]-p2[i];
			d+=delta*delta;
		}
		return sqrt(d);
	}
	double distance_hm(const std::vector<int>& p1, const std::vector<int>& p2)
	{
		double d=0.0;
		for(int i=0;i<p1.size();++i)
			d+= abs(p1[i]-p2[i]);
		return sqrt(d);
	}
private:
	double sigma2_; //parameter for the density kernel
	int nnodes_{0}; //number of items to be clustered
	std::vector<std::vector<double>> distmatrix_; //stores the distance matrix
	std::vector<double> density_; //stores the local density for each item.
	std::vector<int> densityidx_; //stores the indices of the items, sorted by local density, from high to low
	std::vector<std::pair<int,double>> dbetter_; //stores the distance of an item to its closest peer that has a higher local density, and the index of that peer
	std::vector<int> clusterid_; //results of clustering
	std::vector<int> centerindices_;//indices of items at cluster centers
	int nclusters_;
};

}




#endif /* DPCLUSTERING_H_ */
