/*
 * datapaths.cpp
 *
 *  Created on: 2017年4月25日
 *      Author: hyliu
 */

#include "dataio/datapaths.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <cstdlib>
std::string NSPdataio::getenvpath(const std::string & envvar) {
	char *envpath=std::getenv(envvar.c_str());
	std::string res="";
	if(envpath) res=std::string(envpath);
//	char sep= boost::filesystem::path::preferred_separator;
        char sep='/';
	if(res !="" && res.back() != sep) res +=sep;
	return res;
}
std::string NSPdataio::datapath(){
	return getenvpath("ABACUS_DATAPATH");

}
std::string NSPdataio::datafilename(const std::string &filename){
	std::ifstream f(filename.c_str());
	if(f.good()) return filename;
	return datapath()+filename;
}
std::string NSPdataio::downloadedpdbpath(){
	return getenvpath("DOWNLOADED_PDB_PATH");
}
std::string NSPdataio::pdbroot(const std::string &root=""){
	if(root.empty())
		return getenvpath("PDBROOT");
	else
		return root;
}
