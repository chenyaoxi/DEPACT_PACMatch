/*
 * testtopn.cpp
 *
 *  Created on: 2017年1月5日
 *      Author: hyliu
 */

#include "dstl/topn.h"
#include <iostream>
using namespace NSPdstl;

int main() {
	std::vector<double> data{-5,-20,5,8,-17,2,0,-2,-4,100};
	TopN<int> top(5);

	for(int i=0; i<data.size();++i){
		top.push(i,data[i]);
	}
	unsigned int m=top.size();
	for(unsigned int n=0;n<m; ++n) {
		double score;
		std::cout << top.top(&score) ;
		std::cout<<" "<< score<< std::endl;
		top.pop();
	}
}


