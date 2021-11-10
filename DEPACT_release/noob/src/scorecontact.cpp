/*
 * scorecontact.cpp
 *
 *  Created on: 2018年9月4日
 *      Author: hyliu
 */

#include "scorecontact.h"
#include "analyzecontact.h"
#include "atomtypes.h"
#include <cassert>
#include <fstream>
#include <sstream>
#include "dataio/datapaths.h"
using namespace NSPdataio;

using namespace subsitedesign;

std::string copystring(const std::string &s){return s;}

double ContactScoreTable::calcweight(double nobs_tot,double nexp_tot){
	double weight;
	double nscale=nobs_tot>nexp_tot? nobs_tot:nexp_tot;
	static const double NS0{50.0},NS1{5.0};
	if(nscale>=NS0) weight=1.0;
	else { weight=(nscale-NS0)/NS1;
		weight=exp(weight);
	}
	return weight;
}
double ContactScoreTable::getscore(const std::string &ltype, const std::string &ptype, double r,
		double &wght) const{
		std::string lt=ltype;
		std::string pt=ptypecode(ptype);
		if(pt.empty()) pt=ptype;
		if(ltypecode.find(lt) != ltypecode.end()){
			lt=ltypecode.at(lt);
		}
		auto e=table.find(std::make_pair(lt,pt));
		if(e==table.end()) {
			wght=0.0;
			return 0.0;
		}
		wght=e->second.weight;
		int idx=DistBins::distbins().binid(r);
		double rc=DistBins::distbins().bincenter(idx);
		int nidx=idx+1;
		if(r<rc) nidx=idx-1;
		if(nidx<0 || nidx >=DistBins::distbins().num_bins())return e->second.score.at(idx);
		double w1=(r-rc)/(DistBins::distbins().bincenter(nidx)-rc);
		double pav= w1*exp(-(e->second.score.at(nidx)))+(1-w1)*exp(-(e->second.score.at(idx)));
		return -log(pav);
}
bool ContactScoreTable::readentry(std::istream &is){
	std::string line;
	bool firstline=true;
	Entry entry;
	std::string ltype,ptype;
	while(true){
		std::getline(is,line);
		if(!is.good() || line[0]=='&') break;
		std::stringstream sstr(line);
		if(firstline){
			firstline=false;
			sstr >>ltype >>ptype >>entry.nobs_tot >>entry.nexp_tot;
		} else {
			double r,s,no,ne;
			sstr>>r>>s>>no>>ne;
			entry.r.push_back(r);
			double diff=r-DistBins::distbins().bincenter(entry.r.size()-1);
			assert(abs(diff)<1.e-10);
			entry.score.push_back(s);
			entry.nobs.push_back(no);
			entry.nexp.push_back(ne);
		}
	}
	if(entry.r.empty()) return false;
	entry.weight=calcweight(entry.nobs_tot,entry.nexp_tot);
	table[std::make_pair(ltype,ptype)]=entry;
	return true;
}

const ContactScoreTable &ScoreContact::gettable(int level){
	static std::vector<ContactScoreTable> tables(4);
	static bool initialized{false};
	if(!initialized){
		initialized=true;
		std::string line0;
		const std::map<std::string,AtomType> &
			ltypemap=AtomType::getmap();
		for(int i=0;i<4;++i){
			std::string filename=getenvpath("DEPACT_DATAPATH")+"ContactScoreTable"+std::to_string(i)+".dat";
			std::ifstream ifs(filename);
			if (!ifs.is_open())
			{
				std::cout << "Cannot read ContactScoreTable.dat" << std::endl;
				abort();
			}
			std::getline(ifs,line0);
			while(tables[i].readentry(ifs)){ continue;
			}

			for(auto &e:ltypemap){
				if(i==0 || i== 1)tables[i].ltypecode[e.first]= e.second.codename0;
				else tables[i].ltypecode[e.first]=e.second.codename;
			}
			if(i==0 ||i==2) tables[i].ptypecode=copystring;
			else tables[i].ptypecode=subsitedesign::ptype2;
		}
	}
	return tables.at(level);
}
double ScoreContact::repulsivescore(const std::string &ltype,
		const std::string &ptype,double r){
	bool iscarbon=(ltype.at(0)=='C' || ltype.at(0)=='c');
	double radius=2.65;
	if(iscarbon) radius=2.95;
	if(r>radius) return 0.0;
	return exp((radius-r)/0.1);
}
double ScoreContact::score(const std::string & ltype,
		const std::string &ptype,double r){
	if(r<=DistBins::distbins().rmin()|| r>=DistBins::distbins().rmax()) {
		return repulsivescore(ltype,ptype,r);
	}
	double wgt;
	double s;
	for(int i=0;i<4;++i){
		s=gettable(i).getscore(ltype,ptype,r,wgt);
		if(wgt>0.9) break;
	}
	bool ismetal=(ptype=="MG"||ptype=="ZN" || ptype=="FE"||ptype=="MN" ||ptype=="CA"
			||ltype=="MG"||ltype=="ZN" ||ltype=="FE"||ltype=="MN"||ltype=="CA");
	double erep;
	if(ismetal && wgt >0.001)
		erep=repulsivescore(ltype,ptype,r+0.85);
	else
		erep=repulsivescore(ltype,ptype,r);
//	double alpha=0.2; //r-dependent scaling
	double alpha=1.0;
	double rscale=alpha+(1-alpha)*(3.8-r)/1.8;
	if(rscale>1.0) rscale=1.0;
	if(rscale<alpha) rscale=alpha;
//	return (1-wgt)*erep +wgt*s*rscale;
	if(ptype!="HOH"){
		return (1-wgt)*erep+wgt*s*rscale;
	} else {
		return erep+wgt*s*rscale; //some HOH groups in PDB are too close to ligand
	}
}

