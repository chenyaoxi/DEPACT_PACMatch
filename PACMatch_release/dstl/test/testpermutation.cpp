/*
 * testpermutaion.cpp
 *
 *  Created on: 2016年11月30日
 *      Author: hyliu
 */
#include "dstl/permutation.h"
#include <string>
#include <algorithm>
#include <iostream>

using namespace NSPdstl;

int main(int agrc, char ** argv) {
	std::vector<char> origin;
	std::string str(argv[1]);
	for(unsigned int i=0; i<str.size();++i) {
		origin.push_back(str[i]);
	}
	std::sort(origin.begin(),origin.end());
	unsigned int np=Permutation<char>::factorial(origin.size());
	for( unsigned int k=0; k<np; ++k){
		std::vector<char> res=Permutation<char>::getPermutation(origin,k);
		for(char c:res) std::cout <<c;
		std::cout <<std::endl;
	}
}



