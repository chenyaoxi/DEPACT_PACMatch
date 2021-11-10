/*
 * contactstatistics.cpp
 *
 *  Created on: 2018年8月22日
 *      Author: hyliu
 */
#include "atomcontacts.h"
#include "tmpltssas.h"
#include "atomtypessmarts.h"
#include "dataio/inputlines.h"
#include "geometry/calculators.h"
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/filesystem.hpp>
using namespace myobcode;
using namespace subsitedesign;
using namespace OpenBabel;
namespace BFS = boost::filesystem;

int main(int argc, char **argv) {
	OBConversion obconversion;
	obconversion.SetInFormat("sdf");
	std::string molsdf(argv[1]);
	std::ifstream ifs(molsdf);
	obconversion.SetInStream(&ifs);
	std::ofstream ofs((std::string(argv[2])),std::ios_base::app);
	std::vector<std::string> metals{"_MG_","_MN_","_ZN_","_FE_","_K_"};
//	ofs <<" "; //add a delimiter
//	boost::archive::xml_oarchive oa(ofs,1);
	NSPdataio::TextLines lines;
	std::vector<std::string> processed;
	std::string listfile=std::string(argv[3]);
	bool fexist= BFS::exists(listfile) && !BFS::is_empty(listfile);
	if(fexist){
		NSPdataio::TextLines lines;
		lines.init(listfile);
		for(auto &l:lines.lines()) processed.push_back(l);
	}
	bool start=true;
	std::string last_processed;
	if(!processed.empty()){
		last_processed=processed.back();
		start=false;
	}
	std::ofstream ofs_processed;
	ofs_processed.open(listfile.c_str(),std::ios_base::app);
	int pstart=0;
	while(ifs.good()) {
		OBMol mol;
		obconversion.Read(&mol);
		if(!ifs) break;
		int nhvyatoms = mol.NumHvyAtoms();
		std::string moltitle = mol.GetTitle();
		bool doit=true;
		if(!start){
			if(last_processed.find(moltitle) == std::string::npos) continue;
			start=true;
			pstart=processed.size();
		}
		for(int p=pstart; p<processed.size();++p){
			if(processed[p].find(moltitle) != std::string::npos){
				pstart=p+1;
				doit=false;
				break;
			}
		}
		for(auto &m:metals){
			if(moltitle.find(m) != std::string::npos) doit=false;
		}
		if(!doit) continue;
		std::vector<std::vector<std::string>> atypes = findatomtypes(&mol);
		TmpltContacts tc(moltitle);
		if (!tc.good()) continue;
		if(tc.atomcontacts->size()< atypes.size()) continue;
		std::vector<ContactDetails> dtls=collectdetails(tc,atypes);
	//find atoms covalently connected to the protein,ignore them from statistics
		std::set<int> batoms=extbondedligandatoms(dtls,connectivity(&mol));
		ofs<<mol.GetTitle()<<" "<<nhvyatoms <<" "<<atypes.size();
		for(auto &at:atypes){
			if(at.empty()) ofs <<" Un0";
			else ofs <<" "<<at[0];
		}
		ofs<<std::endl;
		ofs <<batoms.size();
		for(auto a:batoms) ofs <<" "<<a;
		ofs <<std::endl;
	//write to file
		for(auto &d:dtls){
			int latm=d.laid;
			if(batoms.find(latm) != batoms.end()) continue;
			ofs<<d.tostring()<<std::endl;
		}
		ofs <<"ENDMOL"<<std::endl;
		ofs.flush();
		ofs_processed<< moltitle <<std::endl;
		ofs_processed.flush();
	} //molecule
}
