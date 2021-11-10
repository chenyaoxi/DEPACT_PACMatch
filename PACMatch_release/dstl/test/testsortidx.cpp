/*
 * testsortidx.cpp
 *
 *  Created on: 2017年3月30日
 *      Author: hyliu
 */

#include "dstl/sortindex.h"
#include <iostream>

int main(){
	std::vector<double> v1{-5,10,8,-3,-4,9,1};
	std::vector<std::string> v2{"-5","10","8","-3","-4","9","1"};
	NSPdstl::sort12(v1,v2);
	for(int i=0;i<v1.size();++i) {
		std::cout <<"("<<v1[i] <<","<<v2[i]<<") ";
	}
	std::cout<<std::endl;
}


