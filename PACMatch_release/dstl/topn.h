/*
 * scoredtype.h
 *
 *  Created on: 2016年12月8日
 *      Author: hyliu
 */

#ifndef DSTL_TOPN_H_
#define DSTL_TOPN_H_
#include <vector>
#include <queue>
namespace NSPdstl {
template <typename T>
struct ScoredType {
		T object;
		double score;
		bool operator<(const ScoredType<T> & st2) const {
			return score <st2.score;
		}
};
template <typename T>
class TopN : public std::priority_queue <ScoredType<T>> {
public:
	unsigned int N;
	TopN():N(0){;}
	void initN(unsigned int n){N=n;}
	TopN(unsigned int n): N(n){;}
	bool push(const T & obj, double score) {
		if(!keep(score)) return false;
		ScoredType<T> st;
		st.score=score;
		st.object=obj;
		if(this->size()>=N) this->pop();
		std::priority_queue <ScoredType<T>>::push(st);
		return true;
	}
	bool push(const T & obj, double score,T *toremove) {
		if(!keep(score)) return false;
		ScoredType<T> st;
		st.score=score;
		st.object=obj;
		if(this->size()>=N) {
			*toremove=std::priority_queue <ScoredType<T>>::top().object;
			this->pop();
		}
		std::priority_queue <ScoredType<T>>::push(st);
		return true;
	}
	bool keep(double score) {
		if (std::priority_queue <ScoredType<T>>::size() >=N)
			if (score >= std::priority_queue <ScoredType<T>>::top().score) return false;
		return true;
	}
	T topobject(double * score) const {
		const ScoredType<T> &st= std::priority_queue <ScoredType<T>>::top();
		*score=st.score;
		return st.object;
	}
	const T & top(double *score) const {
		const ScoredType<T> &st=std::priority_queue <ScoredType<T>>::top();
		*score=st.score;
		return st.object;
	}
	void clear() {
		while (!std::priority_queue <ScoredType<T>>::empty())
			std::priority_queue <ScoredType<T>>::pop();
	}
};
template<typename T>
std::vector<std::pair<T,double>> topN2vector(TopN<T> &topn){
	std::vector<std::pair<T,double>> res;
	res.resize(topn.size());
	for(int i=topn.size()-1;i>=0;--i){
		double score;
		const T & t=topn.top(&score);
		res[i]=std::pair<T,double>(t,score);
		topn.pop();
	}
	for(auto &r:res){
		topn.push(r.first,r.second);
	}
	return res;
}
}


#endif /* DSTL_TOPN_H_ */
