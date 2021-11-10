/*
 * permutation.h
 *
 *  Created on: 2016年11月30日
 *      Author: hyliu
 */

#ifndef DSTL_PERMUTATION_H_
#define DSTL_PERMUTATION_H_
#include <vector>

namespace NSPdstl {

template <typename T>
class Permutation {
public:
	static unsigned int factorial(unsigned int n) {
		{
			unsigned int res = 1;
			for (unsigned int i = 2; i <= n; i++)
				res *= i;
			return res;
		}
	}
	static std::vector<T> getPermutation(const std::vector<T> & original, unsigned int k){
		unsigned int n=original.size();
		std::vector<T> result;
		std::vector<T> reduced = original;
		for(unsigned int i=0; i<n; ++i) {
			result.push_back(helper(reduced,k));
		}
		return result;
	}

private:
	static T helper(std::vector<T> & reduced, unsigned int & k){
		  unsigned int tmp = factorial(reduced.size()-1), i = k/tmp;
		  T res=reduced[i];
		  reduced.erase(reduced.begin()+i);
		  k -= i*tmp;
		  return res;
	}
};
}
#endif /* DSTL_PERMUTATION_H_ */
