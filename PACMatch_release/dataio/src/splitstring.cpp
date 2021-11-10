/*
 * splitstring.cpp
 *
 *  Created on: 2016年11月4日
 *      Author: hyliu
 */

#include "dataio/splitstring.h"
using namespace NSPdataio;
std::vector<int> NSPdataio::integersInString(const std::string & str){
	std::vector<int> substarts=splitString(str,DigitCutter());
	int ns= substarts.size();
	std::vector<int> result;
	substarts.push_back(str.length()+1);
	for(int i=0; i<ns;++i ){
		if(isdigit(str[substarts[i]]))
		result.push_back(std::stoi(str.substr(substarts[i],substarts[i+1]-substarts[i])));
	}
	return result;
}
std::vector<std::string> NSPdataio::alphaWordsInString(const std::string & str){
	std::vector<int> substarts=splitString(str,AlphaCutter());
	int ns= substarts.size();
	std::vector<std::string> result;
	substarts.push_back(str.length()+1);
	for(int i=0; i<ns;++i ){
		if(isalpha(str[substarts[i]]))
		result.push_back(str.substr(substarts[i],substarts[i+1]-substarts[i]));
	}
	return result;
}
std::vector<std::string> NSPdataio::wordsInString(const std::string & str,const std::string & delim){
	DelimCutter cutter(delim);
	std::vector<int> substarts=splitString(str,cutter);
	int ns= substarts.size();
	std::vector<std::string> result;
	substarts.push_back(str.length()+1);
	for(int i=0; i<ns;++i ){
		if(!cutter.isDelim(str[substarts[i]]))
		result.push_back(str.substr(substarts[i],substarts[i+1]-substarts[i]));
	}
	return result;
}



