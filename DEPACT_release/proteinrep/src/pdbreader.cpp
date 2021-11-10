/*
 * pdbreader.cpp
 *
 *  Created on: 2016年11月16日
 *      Author: hyliu
 */

#include "proteinrep/pdbreader.h"
#include "dataio/inputlines.h"
#include "dataio/datapaths.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

using namespace NSPproteinrep;
namespace BFS = boost::filesystem;
void PdbReader::readpdb(const std::string & filename) {
	NSPdataio::TextLines lines;
	lines.init(filename);
	readpdb(lines.lines());
}

void PdbReader::readpdb(std::vector<std::string> & lines) {
	for (auto & line : lines) {
		if (line.substr(0, 6) != "ATOM  " && line.substr(0, 6) != "HETATM")
			continue;
		addRecord(PdbRecord(line));
	}
	mappdbkeyint_=std::shared_ptr<MapPdbKeyInt>(new MapPdbKeyInt(records_));
	for( auto &c:records_) {
		aminoacidsequences_.insert(std::make_pair(c.first,std::vector<std::string>()));
		std::vector<std::string>&seq=aminoacidsequences_.at(c.first);
		int posi=0;
		for(auto &r:c.second){
//			if(r.second[0].label != "ATOM") continue;
			assert(posi++ == mappdbkeyint_->posiNumber(r.first,c.first));
			seq.push_back(r.second[0].residuename);
		}
	}
}

void PdbReader::addRecord(const PdbRecord & record) {
	if (record.failure)
		failure_ = true;
	char chainid = record.chainid;
	if (records_.find(chainid) == records_.end())
		records_.insert(std::make_pair(chainid, ResMapType()));
	ResMapType & resmap = records_.at(chainid);
	ResKeyType reskey = std::make_pair(record.residueid, record.insertionid);
	if (resmap.find(reskey) == resmap.end())
		resmap.insert(std::make_pair(reskey, std::vector<PdbRecord>()));
	resmap.at(reskey).push_back(record);
}

MapPdbKeyInt::MapPdbKeyInt(const typename PdbReader::RecordsType & records) {
	mapchainidint_.init(records);
	mapreskeyint_.resize(records.size(),
			NSPdstl::MapKeyInt<typename PdbReader::ResKeyType>());
	for (auto & c : records)
		mapreskeyint_.at(mapchainidint_.keynumber(c.first)).init(c.second);
}
std::string pdbfilename(const PDBModelID &mid) {
	std::string res;
//	res=boost::algorithm::to_upper_copy(mid.pdbid)+".pdb";
	res="pdb"+mid.pdbid+".ent";
/*	res.resize(mid.pdbid.size());
	std::transform(mid.pdbid.begin(), mid.pdbid.end(), res.begin(), ::tolower);
	res += ".pdb" + std::to_string(mid.biounit);*/
	return res;
}
std::string modelfilename(const PDBModelID &mid) {
	std::string res;
	res = mid.pdbid + "_model" + std::to_string(mid.model) + ".pdb"+std::to_string(mid.biounit);
	return res;
}
static std::vector<std::string>  extractmodellines(const std::string &filename,int modelno){
	std::vector<std::string> res;
	NSPdataio::TextLines lines;
	lines.init(filename);
	bool inmodel = false;
	if (modelno == 0 ||modelno ==1)
		inmodel = true;
	int natoms = 0;
	for (auto &line : lines) {
		if (line.substr(0, 6) != "ATOM  " && line.substr(0, 6) != "HETATM"
				&& line.substr(0, 6) != "ENDMDL"
				&& line.substr(0, 5) != "MODEL") {
			res.push_back(line);
			continue;
		}
		if (line.substr(0, 5) == "MODEL") {
			inmodel = false;
			int num = std::stoi(line.substr(6));
			if (num == modelno)
				inmodel = true;
			continue;
		}
		if (line.substr(0, 6) == "ENDMDL") {
			inmodel = false;
			continue;
		}
		if (inmodel) {
			res.push_back( line);
			++natoms;
		}
	}
	return res;
}
static std::string fndordldpdbfile(const PDBModelID &mid){
	std::string pdbpath = NSPdataio::downloadedpdbpath();
	std::string filename = pdbfilename(mid);
	std::string pfname=pdbpath+filename;
	std::string filenamegz = pdbpath + filename + ".gz";
	BFS::path fpath(pfname);
	BFS::path gzpath(filenamegz);
	bool fexist= BFS::exists(fpath) && !BFS::is_empty(fpath);
	bool gzexist= BFS::exists(gzpath) && !BFS::is_empty(gzpath);
	bool dodld=!fexist && !gzexist;
	if (dodld) {
		std::string cmd;
		cmd="wget -q ftp://ftp.pdbj.org/pub/pdb/data/structures/divided/pdb/"+
				mid.pdbid.substr(1,2)+"/"
//		cmd ="wget -q http://www.rcsb.org/pdb/files/"
				+ filename + ".gz";
/*				"wget -q ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/all/"
						+ filename + ".gz";*/
//#ifndef NDEBUG
		std::cout << "downloading " << filename << ".gz from internet...\n";
//#endif
		FILE *pp = popen(cmd.c_str(), "w");
		if (pp == NULL) {
			std::cout << "downloading unsuccessful\n";
			return std::string();
		}
		std::fclose(pp);
		if (filenamegz != (filename + ".gz")) {
			std::cout <<"moving downloaded file \n";
			cmd = "mv " + filename + ".gz " + filenamegz;
			pp=popen(cmd.c_str(), "w");
			std::fclose(pp);
		}
	}
	if(!fexist){
		std::string cmd;
		cmd = "gunzip " + filenamegz;
		std::cout <<cmd<<std::endl;
		FILE *pp=popen(cmd.c_str(), "w");
		std::fclose(pp);

		// test
		std::cout << filenamegz << " is read" << std::endl;
	}
	if(!BFS::exists(fpath) || BFS::is_empty(fpath)){
		std::cout << "Could not obtain pdbfile "<< filename<<std::endl;
		return std::string();
	}

	// test
	std::cout << "finish " << pfname <<std::endl;

	return pfname;
}
std::string NSPproteinrep::extractpdbmodel(const PDBModelID & mid) {
	std::string filename=fndordldpdbfile(mid);
	return extractpdbmodel(filename, mid);
}

std::string NSPproteinrep::extractpdbmodel(const std::string & filename,
		const PDBModelID &mid) {
	std::vector<std::string> lines=extractmodellines(filename,mid.model);
	if(lines.empty()) return std::string();
	std::string mfile = modelfilename(mid);
	std::ofstream ofs;
	ofs.open(mfile.c_str());
	for(auto &l:lines) ofs <<l<<std::endl;
	return mfile;
	ofs.close();
}
std::shared_ptr<PdbReader> NSPproteinrep::readpdbmodel(const std::string & pdbid,int biounit,
		int modelno){
	PDBModelID mid;
	mid.pdbid=pdbid;
	mid.biounit=biounit;
	mid.model=modelno;
	std::string pdbfile=fndordldpdbfile(mid);
//	if(pdbfile.empty()) return nullptr;
	bool fexist= BFS::exists(pdbfile) && !BFS::is_empty(pdbfile);
	if (!fexist)
	{
		// test
		std::cout << "file " << pdbfile << " not exist" <<std::endl;

		return nullptr;
	}
	std::vector<std::string> lines=extractmodellines(pdbfile,modelno);
	std::string cmd = "gzip " + pdbfile;
	FILE *pp=popen(cmd.c_str(), "w");
	std::fclose(pp);
	std::shared_ptr<PdbReader> reader=std::shared_ptr<PdbReader> (new PdbReader());
	reader->readpdb(lines);
	if (reader->isfailure())
	{
		// test
		std::cout << "reader is failure" << std::endl;

		return nullptr;
	}
	return reader;
}
