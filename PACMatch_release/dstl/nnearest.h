/*
 * nnearest.h
 *
 *  Created on: 2016年3月11日
 *      Author: hyliu
 */

#ifndef NNEAREST_H_
#define NNEAREST_H_
#include <vector>
#include <algorithm>
#include <queue>

namespace domaintree {

template <typename T=long>
class NNearest {
public:
	NNearest(int n=1,double maxd2=1000000000.0): n_(n),maxd2_(maxd2)
	{if (maxd2_ < 0.0) maxd2_=1000000000.0;}


	bool insert(T elem, double dis2, double *maxd2) {
		if(dis2 > maxd2_) return false;
		if( neighbors_.size() < n_) {
			neighbors_.push_back(std::make_pair(elem,dis2));
			if(neighbors_.size() == n_) {
				sort();
				maxd2_=neighbors_[n_-1].second;
				*maxd2=maxd2_;
			}
			return true;
		}
		neighbors_[n_-1]=std::make_pair(elem,dis2);
		sort();
		maxd2_=neighbors_[n_-1].second;
		*maxd2=maxd2_;
		return true;
	}
	double maxd2() {return maxd2_;}
	std::vector<std::pair<T,double>> & neighbors() {return neighbors_;}
	const std::vector<std::pair<T,double>> & neighbors() const {return neighbors_;}
	T nearestone() const {return neighbors_[0].first;}
private:
	int n_;
	std::vector<std::pair<T,double>> neighbors_;
	double maxd2_;
	void sort() {
		std::sort(neighbors_.begin(), neighbors_.end(),
					[](std::pair<T,double> e1,std::pair<T,double> e2)
					->bool {return e1.second < e2.second;});
	}
};
template<typename T>
class QNNearest {
public:
	class Closer {
	public:
		bool operator()(const std::pair<T,double> & a, const std::pair<T,double> &b) {
			return a.second < b.second;
		}
	};
	QNNearest(int n=1,double maxd2=1000000000.0): n_(n),maxd2_(maxd2)
	{if (maxd2_ < 0.0) maxd2_=1000000000.0;}

	bool insert(T elem, double dis2, double *maxd2) {
		if(dis2 > maxd2_) return false;
		if( neighbors_.size() < n_) {
			neighbors_.push(std::make_pair(elem,dis2));
			if(neighbors_.size() == n_) {
				maxd2_=neighbors_.top().second;
				*maxd2=maxd2_;
			}
			return true;
		}
		neighbors_.pop();
		neighbors_.push(std::make_pair(elem,dis2));
		maxd2_=neighbors_.top().second;
		*maxd2=maxd2_;
		return true;
	}
	double maxd2() {return maxd2_;}
	std::vector<std::pair<T,double>> neighbors() {
		std::vector<std::pair<T,double>> nbrs;
		while(!neighbors_.empty()){
			nbrs.push_back(neighbors_.top());
			neighbors_.pop();
		}
		return nbrs;
	}
//	const std::vector<std::pair<T,double>> & neighbors() const {return neighbors_;}
	T nearestone() const {return -1;}
private:
	int n_;
	std::priority_queue<std::pair<T,double>,
	std::vector<std::pair<T,double>>, Closer> neighbors_;
	double maxd2_;
};

}



#endif /* NNEAREST_H_ */
