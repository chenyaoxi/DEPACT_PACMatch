/*
 * scorecontact.h
 *
 *  Created on: 2018年9月4日
 *      Author: hyliu
 */

#ifndef SCORECONTACT_H_
#define SCORECONTACT_H_
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>
#include <cassert>
namespace subsitedesign{

struct ContactScoreTable{
	struct Entry{
		double nobs_tot,nexp_tot;
		std::vector<double> r;
		std::vector<double> nobs;
		std::vector<double> nexp;
		std::vector<double> score;
		double weight;
	};
	std::map<std::pair<std::string,std::string>,Entry> table;
	std::map<std::string,std::string> ltypecode;
	std::string (*ptypecode)(const std::string &pname);
	double getscore(const std::string &ltype, const std::string &ptype, double r,
			double &wght) const;
	bool readentry(std::istream &is);
	static double calcweight(double nobs_tot, double nexp_tot);
};
class ScoreContact{
public:
	static double score(const std::string & ltype, const std::string &ptype,double r);
// scaled score for ligand-mediating water/ion  interactions
	static double score_scaled(const std::string &ltype,const std::string &mtype,double r){
		static std::map<std::string,double> scale{{"HOH",6.5}, {"ZN",5.8},
			{"MG",30.0}, {"MN",11.9}, {"CA",17.8}}; // original HOH:6.0 is too small.
		double s=score(ltype,mtype,r);
		if (scale.find(mtype) == scale.end()) {
			std::cout << "new mediate type: " << mtype << std::endl;
			return 100.0;
		}
		return s/scale[mtype];
	}
// switching value for protein mediating water/ion interactions
	static double theta_m(const std::string &mtype, double mp_score){
		static std::map<std::string,double> EM{{"HOH",-4.0},{"ZN",-4.0},
			{"MG",-5.0},{"MN",-5.5}, {"CA",-8.0}};
		static std::map<std::string,double> SM{{"HOH",3.74},{"ZN",12.2},
			{"MG",13.6},{"MN",11.7}, {"CA",16.0}};
		double x=exp(-(mp_score-EM[mtype])/SM[mtype]);
		return x/(1.0+x);
	}
	static double mp_cutoffs(const std::string &mtype){
		static std::map<std::string,double> cutoffs{{"HOH",-1.0},{"ZN",-1.0},
			{"MG", -1.0},{"MN",-1.0}, {"CA",-1.0}}; // mp_cutoff is for residuepair
		if(cutoffs.find(mtype)== cutoffs.end()) return -1000000.0;
		else return cutoffs[mtype];
	}
private:
	static const ContactScoreTable &gettable(int level);
	static double repulsivescore(const std::string &ltype,
			const std::string &ptype,double r);
};
struct DistBins{
	std::vector<double> rbounds;
	std::vector<double> wghts;
	DistBins(const std::vector<double> &rbds):
		rbounds(rbds){
		wghts.resize(rbounds.size()-1,0.0);
		double wtot=0.0;
		for(int i=0;i<rbounds.size()-1;++i){
			wghts[i]=pow(rbounds[i+1],3)-pow(rbounds[i],3);
			wtot+=wghts[i];
		}
		for(auto &w:wghts) w=w/wtot;
	}
	static const DistBins & distbins(){
//		static DistBins bins({2.1,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.3});
		static DistBins bins({1.5,1.6,1.7,1.8,1.9,2.0,
			2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,3.0,3.2,3.4,3.6,3.8,4.1,4.4,4.7,5.0});
		return bins;
	}
	double rmin() const{return rbounds[0];}
	double rmax() const{return rbounds.back();}
	int num_bins() const {return rbounds.size()-1;}
	int binid(double r) const {
		for(int i=1;i<rbounds.size()-1;++i){
			if(r<rbounds[i]) return i-1;
		}
		return rbounds.size()-2;
	}
	double bincenter(int i) const {
		return (rbounds[i]+rbounds[i+1])*0.5;
	}
};
}




#endif /* SCORECONTACT_H_ */
