/*
 * calcrefdistr.cpp
 *
 *  Created on: 2018年8月23日
 *      Author: hyliu
 */

#include "analyzecontact.h"
#include "atomtypessmarts.h"
#include "scorecontact.h"
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include "openbabel/math/vector3.h"
#include "geometry/xyz.h"
using namespace OpenBabel;
using namespace myobcode;
using namespace subsitedesign;
int main(int argc,char **argv){
	OBConversion obconversion;
	obconversion.SetInFormat("sdf");
	std::string molsdf(argv[1]);
	std::ifstream ifs(molsdf);
	obconversion.SetInStream(&ifs);
	std::map<std::string,std::vector<double>> distr1;
	std::map<std::string,std::vector<double>> distr2;
	std::map<std::string,double> atype_count;
	int nmol=0;
	int nmolmax=std::atoi(argv[2]);
	int nskipmol=std::atoi(argv[3]);
	int nskipped=0;
	while(ifs.good()) {
			if(nmol >= nmolmax) break;
			OBMol mol;
			obconversion.Read(&mol);
			if(!ifs) break;
			int nhvyatoms = mol.NumHvyAtoms();
			if(nhvyatoms<5) continue;
			nskipped++;
			if(nskipped<=nskipmol) continue;
			auto &amap=AtomType::getmap();
			std::vector<std::vector<std::string>> atypes = findatomtypes(&mol);
			int nnco=0;
			std::vector<std::string> at_simp;
			for(auto&a:atypes){
				if(a.empty()) at_simp.push_back("Un0");
				else at_simp.push_back(amap.at(a[0]).codename);
			}
			for(auto &a:at_simp){
				if(a[0]=='H') continue;
				if(atype_count.find(a)== atype_count.end())
					atype_count[a]=0.0;
				atype_count[a] +=1.0;
			}
			std::vector<NSPgeometry::XYZ> crds;
			for (int i = 0; i < mol.NumAtoms(); ++i) {
				vector3 v3=mol.GetAtom(i+1)->GetVector();
				crds.push_back(NSPgeometry::XYZ(v3.x(),v3.y(),v3.z()));
			}
			samplerefdist(at_simp, crds,
					distr1,
					distr2);
			nmol++;
	}
	double all=0.0;
	for(auto &ac:atype_count){
		all+=ac.second;
	}
	atype_count["all"]=all;
	double ncon_all1=0.0;
	double ncon_all2=0.0;
	for(auto d:distr1["all"]) ncon_all1+=d;
	for(auto d:distr2["all"]) ncon_all2+=d;
	for(int i=0;i<2;i++){
		std::map<std::string,std::vector<double>> *distr=&distr1;
		double ncon_all=ncon_all1;
		if(i==1) {
			distr=&distr2;
			ncon_all=ncon_all2;
		}
		for(auto &md:*distr){
			std::cout <<"&"<<std::endl;
			std::cout <<"#"<<md.first;
			double ncon_tot=0.0;
			for(auto d:md.second) ncon_tot+=d;
			double ncon_ratio=(atype_count[md.first]/atype_count["all"])
					*ncon_all/ncon_tot;
			std::cout <<" "<<1.0/ncon_ratio<<std::endl;
			for(int bidx=0;bidx<md.second.size();++bidx){
				std::cout <<DistBins::distbins().bincenter(bidx)<<" "
						<<md.second[bidx]/ncon_tot<<std::endl;
			}
		}
	}
}


