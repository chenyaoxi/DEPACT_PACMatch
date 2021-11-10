/*
 * pdbfolder.cpp
 *
 *  Created on: 2016年2月23日
 *      Author: hyliu
 */

#include <stdlib.h>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>   // includes all needed Boost.Filesystem declarations
#include <iostream>               // for std::cout
#include "dataio/pdbfolder.h"
#define FOLDERSEP "/"
namespace BFS=boost::filesystem;
using namespace NSPdataio;
std::string PdbFolder::subfolder(std::string name) {
	std::string dirname = root_ +FOLDERSEP+ name +FOLDERSEP;
    BFS::path pdbdir(dirname);
    if(!BFS::exists(pdbdir)) {
    	if(BFS::create_directories(pdbdir))
    		std::cout <<"Cannot make directory " << dirname <<"\n";
    		return std::string();
    	}
    return dirname;
}

std::string PdbFolder::getfile(std::string pdbid, const std::string & ext){
		boost::algorithm::to_lower(pdbid);
		std::string folder=subfolder(pdbid.substr(2,1));
	if(folder.empty()) return std::string();
	std::string pdbfilename = folder +"pdb"+ pdbid + "." + ext;
	BFS::path pdbfile(pdbfilename);
//	std::cout <<pdbfilename <<std::endl;
	if(BFS::exists(pdbfile) && !BFS::is_empty(pdbfile)) return pdbfilename;
	if(ext != "ent") {
		std::cout <<" file " <<pdbfilename << " none exist or empty\n";
		return std::string();
	}
	return downloadpdb(pdbid, pdbfilename);
}

std::string NSPdataio::downloadpdb(std::string pdbid,std::string localfile) {
	std::string cmd;
	cmd = "wget -q -O - http://www.rcsb.org/pdb/files/"
			+ boost::algorithm::to_upper_copy(pdbid)+".pdb.gz"
			+" | gunzip -c - > " + localfile;
#ifndef NDEBUG
	std::cout<<"retrieving pdb " <<pdbid <<"from internet...\n";
#endif
	FILE *pp=popen(cmd.c_str(),"w");
	if(pp == NULL) {
		std::cout <<"get pdb unsuccessful\n";
		return std::string();
	}
	pclose(pp);
	BFS::path file(localfile);
	if(BFS::exists(file) && !BFS::is_empty(file)) {
#ifndef NDEBUG
		std::cout <<"pdbfile retrieved as " <<localfile <<"\n";
#endif
		return localfile;
	} else {
		std::cout <<"get pdb file unsuccessful\n";
		return std::string();
	}
}


