/*
 * testsplitstring.cpp
 *
 *  Created on: 2016年11月4日
 *      Author: hyliu
 */

#include "dataio/splitstring.h"
#include <iostream>

using namespace NSPdataio;
int main(){
	std::string teststr{"This 56789test,random;     567 ,;;t1"};
	std::cout <<"TestString="<<teststr <<std::endl;

	std::vector<int> ints=integersInString(teststr);
	for(int i:ints) std::cout <<i <<"\t";
	std::cout <<std::endl;

	std::vector<std::string> awords=alphaWordsInString(teststr);
	for(auto & w:awords) std::cout <<w <<"\t";
	std::cout <<std::endl;


	std::vector<std::string> words=wordsInString(teststr,",;");
	for(auto & w:words) std::cout <<w <<"\t";
	std::cout <<std::endl;
}


