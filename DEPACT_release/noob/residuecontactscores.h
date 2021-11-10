/*
 * tmpltresiduescores.h
 *
 *  Created on: 2018年9月4日
 *      Author: hyliu
 */

#ifndef RESIDUECONTACTSCORES_H_
#define RESIDUECONTACTSCORES_H_
#include "atomcontacts.h"

namespace subsitedesign{

//Summarized contact scores
struct ResidueContactScores{
	typedef std::pair<int,int> ResidueID;
	//scores of each atom in each protein residue versus each ligand atom
	std::map<ResidueID,std::map<std::string,std::vector<double>>> plscores;
	//scores of mediating atom versus each atom in each protein residue
	std::map<ResidueID,std::map<ResidueID,std::map<std::string,double>>> mpscores;
	//scores of mediating atom versus each each ligand atom
	std::map<ResidueID,std::vector<double>> mlscores;
	ResidueContactScores(){;};
	void insertplscore(ResidueID prid,std::string paname,int nlatom,
			int laidx,double s){
		if(plscores.find(prid)==plscores.end()){
			plscores[prid]=std::map<std::string,std::vector<double>>();
		}
		if(plscores[prid].find(paname)==plscores[prid].end()){
			plscores[prid][paname]=std::vector<double>(nlatom,0.0);
		}
		plscores[prid][paname][laidx]=s;
	}
	void insertmlscore(ResidueID mrid,int nlatom,
			int laidx,double s){
		if(mlscores.find(mrid)==mlscores.end()){
			mlscores[mrid]=std::vector<double>(nlatom,0.0);
		}
		mlscores[mrid][laidx]=s;
	}
	void insertmpscore(ResidueID mrid,ResidueID prid,
			std::string paname,double s){
		if(mpscores.find(mrid)==mpscores.end()){
			mpscores[mrid]=std::map<ResidueID,std::map<std::string,double>>();
		}
		if(mpscores[mrid].find(prid)==mpscores[mrid].end()){
			mpscores[mrid][prid]=std::map<std::string,double>();
		}
		mpscores[mrid][prid][paname]=s;
	}
	ResidueContactScores(const TmpltContacts &tcons,
			const std::vector<std::string> &atypes);
	void setup(const std::vector<AtomContacts> &acons,
			const std::vector<std::vector<NSPproteinrep::AAConformer>> &aas,
			const std::vector<std::string> &latypes);
	double totalplscore(ResidueID r) const {
		if(plscores.find(r)== plscores.end()) return 0.0;
		double sum=0.0;
		for(auto & ms:plscores.at(r)) { ///<every atom in residue
			for(auto s:ms.second) ///<every ligand atom
				sum+=s;
		}
		return sum;
	}
	double backboneplscore(ResidueID r) const {
		static std::set<std::string> backboneatoms{"N","CA","C","O"};
		if(plscores.find(r)== plscores.end()) return 0.0;
		double sum=0.0;
		for(auto & ms:plscores.at(r)) { ///<every atom in residue
			if(backboneatoms.find(ms.first)!=backboneatoms.end())
			for(auto s:ms.second) ///<every ligand atom
				sum+=s;
		}
		return sum;
	}
	double sidechainplscore(ResidueID r) const {
		static std::set<std::string> backboneatoms{"N","CA","C","O"};
		if(plscores.find(r)== plscores.end()) return 0.0;
		double sum=0.0;
		for(auto & ms:plscores.at(r)) { ///<every atom in residue
			if(backboneatoms.find(ms.first)==backboneatoms.end())
			for(auto s:ms.second) ///<every ligand atom
				sum+=s;
		}
		return sum;
	}
	double totalplscore(ResidueID r, std::string a) const{
		if(plscores.find(r)== plscores.end()) return 0.0;
		if(plscores.at(r).find(a)== plscores.at(r).end()) return 0.0;
		double sum=0.0;
		for(auto s:plscores.at(r).at(a)) sum+=s; ///<every ligand atom
		return sum;
	}
	double mpscore(ResidueID m,ResidueID r)const {
		if(mpscores.find(m)==mpscores.end()) return 0.0;
		if(mpscores.at(m).find(r)==mpscores.at(m).end()) return 0.0;
		double sum=0.0;
		for(auto s:mpscores.at(m).at(r)) sum+=s.second;
		return sum;
	}
	double backbonempscore(ResidueID m,ResidueID r)const {
		static std::set<std::string> backboneatoms{"N","CA","C","O"};
		if(mpscores.find(m)==mpscores.end()) return 0.0;
		if(mpscores.at(m).find(r)==mpscores.at(m).end()) return 0.0;

		double sum=0.0;
		for(auto s:mpscores.at(m).at(r)) {
			if(backboneatoms.find(s.first) != backboneatoms.end())
				sum+=s.second;
		}
		return sum;
	}
	double sidechainmpscore(ResidueID m,ResidueID r)const {
		static std::set<std::string> backboneatoms{"N","CA","C","O"};
		if(mpscores.find(m)==mpscores.end()) return 0.0;
		if(mpscores.at(m).find(r)==mpscores.at(m).end()) return 0.0;

		double sum=0.0;
		for(auto s:mpscores.at(m).at(r)) {
			if(backboneatoms.find(s.first) == backboneatoms.end())
				sum+=s.second;
		}
		return sum;
	}
	double totalmpscore(ResidueID m) const {
		if(mpscores.find(m)==mpscores.end()) return 0.0;
		double sum=0.0;
		for(auto &ps:mpscores.at(m)) { ///<every  protein residue
			for(auto s:ps.second) sum +=s.second;  ///<every atom in protein residue
		}
		return sum;
	}
	double totalmlscore(ResidueID m) const {
		if(mlscores.find(m)==mlscores.end()) return 0.0;
		double sum=0.0;
		for(auto s:mlscores.at(m))sum+=s;  ///<every ligand atom
		return sum;
	}
	void print(std::ostream &os) const {
		for(auto &pl:plscores){
			os<<"plscore "<<pl.first.first<<"-"<<pl.first.second<<": "<<
					totalplscore(pl.first)<<std::endl;
		}
		for(auto &ml:mlscores){
			os<<"mlscore "<<ml.first.first<<"-"<<ml.first.second<<": "<<
					totalmlscore(ml.first) << std::endl;
		}
		for(auto &ml:mpscores){
			os<<"mpscore "<<ml.first.first<<"-"<< ml.first.second<<": "<<
					totalmpscore(ml.first)<< std::endl;
		}
	}
};

}



#endif /* RESIDUECONTACTSCORES_H_ */
