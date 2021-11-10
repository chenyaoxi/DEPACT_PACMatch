/*
 * analyzecontacts.cpp
 *
 *  Created on: 2018年8月22日
 *      Author: hyliu
 */
#include "analyzecontact.h"
#include "atomcontacts.h"
#include "atomtypessmarts.h"
#include "scorecontact.h"
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <sstream>
using namespace myobcode;
using namespace subsitedesign;
/*
 * reasoning of statistical scores
 * expect_contact_freq=total_contact_count/N_ltotal*N_ptotal;
 * total_pair_contact_expected=expect_contact_freq*(N_latomtype*N_ptype)
 * pair expected_in_bin= w_bin*total_contact_count_expected
 * w_bin=(V_binshell)/sum_bin(V_binshell)
 * V_binshell=4/3*pi*(Rmax_bin^3-Rmin_bin^3)
 * ene_pair_bin=-log(pair_observed_in_bin/pair_expected_in_bin)
 * Score ligand-protein, ligand-mediator interactions:
 * 		ligand->protain+mediating
 * 		N_ltotal=N-ligand; N_ptotal=N_protein+N_mediating
 * Score mediator-protein interactions
 *    ligand+mediating->protein
 *    N_ltotal=N-ligand+N_mediating N_ptotal=N_protein
 *
 */

struct RefDistr{
	std::vector<double> wghts;
	double rate;
	double readtimes;
	RefDistr():wghts(DistBins::distbins().num_bins(),0.0),rate(0.0),readtimes(0.0){;}
	void average(){
		for(auto &w:wghts) w =w/readtimes;
		rate=rate/readtimes;
		readtimes=1.0;
	}
	void print(std::ostream &ofs){
		ofs <<" "<<rate <<std::endl;
		for(int i=0;i<wghts.size();++i){
			ofs<<DistBins::distbins().bincenter(i)<<" "<<wghts[i]<<std::endl;
		}
	}
};
void readdistrs(std::istream &ifs,
		std::map<std::string,RefDistr> & distrs){
	std::string line;
	bool indist=false;
	std::string atype;
	double ratio;
	while(ifs.good()){
		std::getline(ifs,line);
		if(!ifs.good()) break;
		if(line[0]=='&'){
			std::getline(ifs,line);
			std::string st=line.substr(1);
			std::stringstream sstr(st);
			sstr>>atype;
			sstr>>ratio;
			if(distrs.find(atype)==distrs.end()){
				distrs.insert(std::make_pair(atype,RefDistr()));
			}
			distrs[atype].rate +=ratio;
			distrs[atype].readtimes+=1.0;
			continue;
		}
		std::stringstream sstr(line);
		double r,w;
		sstr >> r >>w;
		int bin=DistBins::distbins().binid(r);
		distrs[atype].wghts[bin] +=w;
	}
}

struct MolContactSum {
	typedef std::pair<int, int> ResID;
	//all contacts made by a protein atom(inner key) in a protein residue (outter key)
	typedef std::map<ResID,
			std::map<std::string, std::vector<const ContactDetails *>>>SetSum;
	SetSum dcsetsum; //summary of direct contacts
	std::map<ResID,SetSum> indcsetsums;//summary indirect contacts, ResId is that of mediator
	SetSum allsetsums;
	std::map<ResID,std::map<int,double>> inddistances;//mediator -ligand atom distances
	std::shared_ptr<std::vector<ContactDetails>> directdtls;
	std::shared_ptr<std::map<ResID,std::vector<ContactDetails>>>indirectdtls;
	MolContactSum(const std::vector<ContactDetails> & dtls) {
		directdtls=std::shared_ptr<std::vector<ContactDetails>> (new std::vector<ContactDetails>());
		indirectdtls=std::shared_ptr<std::map<ResID,std::vector<ContactDetails>>>(
				new std::map<ResID,std::vector<ContactDetails>>());
		for(auto &d:dtls) {
			if(d.indirect==0) {
				directdtls->push_back(d);
			} else {
				ResID mr(d.mcid,d.mrid);
				if(indirectdtls->find(mr)==indirectdtls->end())
				(*indirectdtls)[mr]=std::vector<ContactDetails>();
				(*indirectdtls)[mr].push_back(d);
				if(inddistances.find(mr)==inddistances.end())
				inddistances[mr]=std::map<int,double>();
				inddistances[mr][d.laid]=d.mdistance;
			}
		}
		allsetsums=anamol(dtls);
		dcsetsum=anamol(*directdtls);
		for(auto &md:*indirectdtls) {
			if(indcsetsums.find(md.first)==indcsetsums.end())
			indcsetsums[md.first]=anamol(md.second);
		}
	}
	static SetSum anamol(const std::vector<ContactDetails> &dts) {
		SetSum sum;
		for (auto &c:dts) {
			std::pair<int,int> residue(c.pcid,c.prid);
			if(sum.find(residue) == sum.end()) {
				sum[residue]=std::map<std::string,std::vector<const ContactDetails*>>();
			}
			auto &entryr=sum[residue];
			if(entryr.find(c.paname)== entryr.end()) {
				entryr[c.paname]=std::vector<const ContactDetails*>();
			}
			entryr[c.paname].push_back(&c);
		}
		return sum;
	}
};

double calcscore(double total_expect,double expect,double observe){
	double w=20/total_expect;
	w=w*w*w;
	double wmin;
	if(observe>=1 && observe/expect >=5) {
		observe=observe-0.8;
	}
	return (observe+w*expect)/((1.0+w)*expect);
}
int main(int argc, char** argv) {

	std::map<std::string, double> num_ligandatoms;
	std::map<std::string, double> num_platoms;
	std::map<std::string, std::map<std::string,double>> num_pmatoms;
	std::map<std::string, double> num_mlatoms;

	std::map<std::string, std::vector<double>> num_lpcontacts;
	std::map<std::string, std::vector<double>> num_lmcontacts;
	std::map<std::string, std::map<std::string,std::vector<double>>> num_pmcontacts;
	const DistBins &bins = DistBins::distbins();
	std::ifstream ifs((std::string(argv[1])));
	std::map<std::string, AtomType> map = AtomType::getmap();
	std::string (*getpatype)(const std::string &, std::string);
	std::string (*getlatype)(const std::string &);
	std::map<std::string,std::string> nmmap;
	int mode1=std::atoi(argv[3]);
	int mode2=std::atoi(argv[4]);
	if(mode1!=2) getlatype=subsitedesign::getcodename0;
	else getlatype=subsitedesign::getcodename;
	if(mode2!=2) getpatype=subsitedesign::proteinatomtype;
	else getpatype=subsitedesign::proteinatomtype2;
	while (ifs.good()) {
		int nlatoms;
		std::vector<std::string> atypes;
		std::set<int> excluded;
		std::string molname;
		std::vector<ContactDetails> dts = getdtlsnextmol(ifs, molname, nlatoms, atypes,excluded);
		if (dts.empty() || nlatoms == 0)
			continue;
		int idx=-1;
		for (auto &at : atypes) {
			++idx;
			std::string latype = getlatype(at);
			nmmap[latype]=subsitedesign::getcodename(at);
			if (latype[0] == 'H')
				continue;
			if(excluded.find(idx) != excluded.end()) continue;
			if (num_ligandatoms.find(latype) == num_ligandatoms.end())
				num_ligandatoms[latype] = 0.0;
			num_ligandatoms[latype] += 1.0;
		}
		MolContactSum molsum(dts);
		for (auto &e : molsum.dcsetsum) {
			std::string prname = e.second.begin()->second[0]->prname;
			for (auto &sd : e.second) {
				std::string patname = getpatype(prname, sd.first);
				if(patname.empty()) continue;
				if (num_platoms.find(patname) == num_platoms.end()) {
					num_platoms[patname] = 0.0;
				}
				num_platoms[patname] += 1.0;
				for (auto d : sd.second) {
					std::string latype;
					try{
						latype= getlatype(d->latype);
					} catch (std::exception &e){
						std::cout <<d->latype<<std::endl;
						exit(1);
					}
					std::string con_type = latype + ":" + patname;
					if (num_lpcontacts.find(con_type) == num_lpcontacts.end())
						num_lpcontacts[con_type] = std::vector<double>(
								bins.num_bins(), 0.0);
					int idx = bins.binid(d->distance);
					num_lpcontacts[con_type][idx] += 1.0;
				}
			}
		} //direct contacts
		for (auto &me : molsum.indcsetsums) {
			std::string mname =
					me.second.begin()->second.begin()->second[0]->mligand;
			if (num_mlatoms.find(mname) == num_mlatoms.end()) {
				num_mlatoms[mname] = 0.0;
			}
			num_mlatoms[mname] += 1.0;
			for (auto &id : molsum.inddistances[me.first]) {
/*				if(id.first >=atypes.size()){
					break;
				}
				if(map.find(atypes[id.first])==map.end()) {
					break;
				}*/
				std::string latype;
				try {
					latype=getlatype(atypes[id.first]);
				} catch (std::exception &e){
					std::cout << atypes.size()<< " "<<id.first<<std::endl;
					for(auto &t:atypes) std::cout <<" "<<t;
					std::cout<<std::endl;
					exit(1);
				}
				std::string con_type = latype + ":" + mname;
				if (num_lmcontacts.find(con_type) == num_lmcontacts.end()) {
					num_lmcontacts[con_type] = std::vector<double>(
							bins.num_bins(), 0.0);
				}
				num_lmcontacts[con_type][bins.binid(id.second)] += 1.0;
			}
			if(num_pmatoms.find(mname)==num_pmatoms.end()){
				num_pmatoms[mname]=std::map<std::string,double>();
			}
			if(num_pmcontacts.find(mname)==num_pmcontacts.end()){
				num_pmcontacts[mname]=std::map<std::string,std::vector<double>>();
			}
			for (auto &e : me.second) {
				std::string prname = e.second.begin()->second[0]->prname;
				for (auto &sd : e.second) {
					std::string patname = getpatype(prname, sd.first);
					if(patname.empty()) continue;
					if (num_pmatoms[mname].find(patname) == num_pmatoms[mname].end()) {
						num_pmatoms[mname][patname] = 0.0;
					}
					num_pmatoms[mname][patname] += 1.0;
					for (auto d : sd.second) {
						int idx = bins.binid(d->distance);
						std::string ct = patname;
						if (num_pmcontacts[mname].find(ct) == num_pmcontacts[mname].end()) {
							num_pmcontacts[mname][ct] = std::vector<double>(
									bins.num_bins(), 0.0);
						}
						num_pmcontacts[mname][ct][idx] += 1.0;
					}
				}
			}
		}//indirect
	} //ifs good
	double ntot_latoms=0;
	for(auto &n:num_ligandatoms) ntot_latoms+=n.second;
	double ntot_platoms=0;
	for(auto &n:num_platoms) ntot_platoms+=n.second;
	double ntot_mlatoms=0;
	for(auto &n:num_mlatoms) ntot_mlatoms+=n.second;

	double ntot_lpcontacts=0;
	for(auto &e:num_lpcontacts){
		for(auto &d:e.second) ntot_lpcontacts+=d;
	}
	double ntot_lmcontacts=0;
	for(auto &e:num_lmcontacts){
		for(auto &d:e.second) ntot_lmcontacts+=d;
	}
	double ntot_lcontacts=ntot_lpcontacts+ntot_lmcontacts;
	double ntot_pmlatoms=ntot_platoms+ntot_mlatoms;
	double rate_av= ntot_lcontacts/ntot_latoms;
//	std::map<std::string,double> refratio;
//	std::map<std::string,std::vector<double>> refwghts;
	//for test
/*	for(auto & e:num_ligandatoms){
		refratio[e.first]=1.0;
		refwghts[e.first]=std::vector<double>(bins.num_bins(),1.0);
	}*/

	std::ifstream ifsd((std::string(argv[2])));
	std::map<std::string,RefDistr> distrs;
	readdistrs(ifsd,distrs);
	for(auto &d:distrs){
		d.second.average();
/*		std::cout <<"&"<<std::endl;
		std::cout <<"#"<<d.first;
		d.second.print(std::cout);*/
	}
//	exit(0);
	std::string mtype_equiv="O.31";
	for(auto & e:num_lpcontacts){
		std::string latype=e.first.substr(0,e.first.find(':'));
		std::string patype=e.first.substr(e.first.find(':')+1);
		double nc=0;
		for(auto n:e.second) {
			nc+=n;
		}
		double	nc_expect =rate_av*num_ligandatoms[latype]*
					num_platoms[patype]/ntot_pmlatoms*distrs[nmmap[latype]].rate;
		std::vector<double> score(bins.num_bins(),0);
		std::cout <<"&"<<std::endl;
		std::cout <<latype <<" "<<patype<< " "<<nc<<" "<<nc_expect<<std::endl;
		for(int i=0;i<bins.num_bins();++i){
			double expect=nc_expect*
					distrs[nmmap[latype]].wghts[i]+1.e-5;
			score[i]=calcscore(0.5*(nc_expect+nc),expect,e.second[i]);
			std::cout <<bins.bincenter(i) <<" "<<score[i]
					<<" "<<e.second[i]<<" "<<expect<<std::endl;
		}
	}
	for(auto & e:num_lmcontacts){
		std::string latype=e.first.substr(0,e.first.find(':'));
		std::string matype=e.first.substr(e.first.find(':')+1);
		double nc=0;
		for(auto n:e.second) {
			nc+=n;
		}
		double nc_expect =rate_av*num_ligandatoms[latype]*
					num_mlatoms[matype]/ntot_pmlatoms*distrs[nmmap[latype]].rate;
		std::vector<double> score(bins.num_bins(),0);

		std::cout <<"&"<<std::endl;
		std::cout <<latype <<" "<<matype<<" "<<nc<<" "<<nc_expect<<std::endl;
		for(int i=0;i<bins.num_bins();++i){
			double expect=nc_expect*
					distrs[nmmap[latype]].wghts[i]+1.e-5;
			score[i]=calcscore(0.5*(nc_expect+nc),expect,e.second[i]);
			std::cout <<bins.bincenter(i) <<" "<<score[i]
					<<" "<<e.second[i]<<" "<<expect<<std::endl;
		}
	}
	std::map<std::string,double> ref_pref;
	for(auto & npl:num_platoms){
		ref_pref[npl.first]=npl.second/ntot_platoms;   //protein atomtype distr, ligand
	}
	for(auto &me:num_pmcontacts){
		std::string matype=me.first;
		double ntot_pm=0.0;
		for(auto &n:me.second) {
			for(auto d:n.second) ntot_pm+=d;
		}
		double ntot_pmatoms=0.0;
		for(auto &n:num_pmatoms[me.first]){
			ntot_pmatoms +=n.second;
		}
		for(auto &pe:me.second){
			std::string patype=pe.first;
			std::vector<double> score(bins.num_bins(),0.0);
			double pref=num_pmatoms[matype][patype]/ntot_pmatoms;
			double ncon=0.0;
			for (auto d:pe.second)ncon+=d;
			std::cout<<"&"<<std::endl;
			std::cout<<matype<<" "<<patype<<" "<<ncon<<" "<<ncon<<std::endl;
			for(int i=0;i<bins.num_bins();++i){
				double expect=ncon*distrs[mtype_equiv].wghts[i]+1.e-5;
				score[i]=pref/ref_pref[patype]*calcscore(ncon,expect,pe.second[i]);
				std::cout<<bins.bincenter(i)<<" "<<score[i]
						<<" "<<pe.second[i]<<" "<<expect<<std::endl;
			}
		}
	}
}
