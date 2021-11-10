/*
 * tmpltssas.cpp
 *
 *  Created on: 2018年8月10日
 *      Author: hyliu
 */
#include "tmpltssas.h"
using namespace subsitedesign;
std::vector<std::string> TmpltSSAs::parsetmpltname(const std::string &name){
	std::vector<std::string> res;
	std::vector<int> delimposi;
	delimposi.push_back(-1);
	for(int i=1;i<name.size();++i){
		if(name[i]=='_') delimposi.push_back(i);
	}
	for(int i=0;i<delimposi.size()-1;++i){
		res.push_back(name.substr(delimposi[i]+1,delimposi[i+1]-(delimposi[i]+1)));
	}
	return res;
}
static inline std::string & trim(std::string &str){
	str.erase(str.begin(), std::find_if(str.begin(), str.end(), [](int ch) {
	        return !std::isspace(ch);}));
	str.erase(std::find_if(str.rbegin(), str.rend(), [](int ch) {
	        return !std::isspace(ch);
	    }).base(), str.end());
	return str;
}
TargetStruct::TargetStruct(std::string name,
		const std::vector<std::string> & atypes,std::istream &ispdb){
	molname=name;
	atomtypes=atypes;
	std::string line;
	std::vector<NSPproteinrep::PdbRecord> records;
	while(getline(ispdb,line)){
		if(line.substr(0,6) !="ATOM  " && line.substr(0,6) !="HETATM") continue;
		NSPproteinrep::PdbRecord nrecord(line);
		std::string e;
		e.push_back(nrecord.elementname[0]);
		e.push_back(nrecord.elementname[1]);
		nrecord.atomname=trim(e)
				+std::to_string(nrecord.atomid);
		records.push_back(nrecord);
	}
	assert(records.size()==atomtypes.size());
	conformer=NSPproteinrep::make_aaconformer(records);
}
