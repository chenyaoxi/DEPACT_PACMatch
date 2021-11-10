/*
 * testcontrolfile.cpp
 *
 *  Created on: 2017年8月19日
 *      Author: hyliu
 */
#include "dataio/controlfile.h"

using namespace NSPdataio;

int main(int argc, char ** argv){
	ControlFile &cf1=ControlFileInstances::getstaticinstance();
	cf1.readfile(std::string(argv[1]));
	ControlFile &cf2=ControlFileInstances::getnamedinstance("name1");
	cf2.readfile(std::string(argv[2]));
	auto &cf3=ControlFileInstances::copystaticinstance("name2");
	cf3.readfile(std::string(argv[3]));
	std::cout <<"CF1: --"<<std::endl;
	ControlFileInstances::getinstance().write(std::cout);
	std::cout <<"CF2: --" <<std::endl;
	ControlFileInstances::getnamedinstance("name1").write(std::cout);
	std::cout <<"CF3: --" <<std::endl;
	ControlFileInstances::getinstance("name2").write(std::cout);
}



