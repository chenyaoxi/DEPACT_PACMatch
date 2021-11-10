/*
 * controlfile.cpp
 *
 *  Created on: 2017年8月19日
 *      Author: hyliu
 */
#include "dataio/inputlines.h"
#include "dataio/controlfile.h"
#include <sstream>
using namespace NSPdataio;
using namespace NSPdstl;

void ControlFile::write(std::ostream &os) const{
	for(auto &m:controllines_){
		os<< "START " <<m.first<<std::endl;
		for (auto &line:m.second){
			os << "   "<<line<<std::endl;
		}
		os<< "END   " <<m.first<<std::endl<<"---------------------------------------"<<std::endl;
	}
}
void ControlFile::readfile(std::istream &is){
	char buffer[600];
	bool inblock{false};
	std::string controlname;
	while(true){
		is.getline(buffer,600);
		if(!is.good()) break;
		std::string line(buffer);
		std::stringstream sstr(line);
		std::vector<std::string> words=parseline(line,std::vector<int>());
		if(words.size()==0) continue;
		if(!inblock) {
			if(words[0] == "START" ) {
				controlname=words[1];
				if(definedcontrolnames_.find(controlname)==definedcontrolnames_.end())
					status_.set(readundefined);
				inblock=true;
				if(controllines_.find(controlname) == controllines_.end())
					controllines_.insert(std::make_pair(controlname,std::vector<std::string>()));
			} else if(words[0] == "INCLUDE") {
				readfile(words[1]);
				if(openfilefailed()){
					std::cout <<"Open included file " << words[1] << "failed." <<std::endl;
					abort();
				}
			} else {
				std::cout <<"Ignored input line: " <<line <<std::endl;
			}
		} else {
			if(words[0] == "END"){
				if(controlname != words[1]) {
					std::cout <<"Unmatched names in START and END lines in input control: ";
					std::cout <<controlname <<" v.s. " <<words[1]<<std::endl;
				}
				inblock=false;
				controlname="";
			} else {
				if(line.find('=') == std::string::npos){
					std::cout <<"Invalid control line: " <<line <<std::endl;
					abort();
				}
				controllines_.at(controlname).push_back(line);
			}
		}

	}
	if(inblock) {
		std::cout <<"Last block not ended with an \"END controlname\" line in control input."<<std::endl;
		abort();
	}
}
