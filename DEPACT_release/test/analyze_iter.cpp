/*
 * analyze_iter.cpp
 *
 *  Created on: 2018年9月8日
 *      Author: hyliu
 */

#include "analyzecontact.h"
#include "atomcontacts.h"
#include "atomtypessmarts.h"
#include "scorecontact.h"
#include "dstl/randomengine.h"
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
std::string (*getpatype)(const std::string &, std::string);
std::string (*getlatype)(const std::string &);
struct RefDistr {
	std::vector<double> wghts;
	double rate;
	double readtimes;
	RefDistr() :
			wghts(DistBins::distbins().num_bins(), 0.0), rate(0.0), readtimes(
					0.0) {
		;
	}
	void average() {
		for (auto &w : wghts)
			w = w / readtimes;
		rate = rate / readtimes;
		readtimes = 1.0;
	}
	void print(std::ostream &ofs) {
		ofs << " " << rate << std::endl;
		for (int i = 0; i < wghts.size(); ++i) {
			ofs << DistBins::distbins().bincenter(i) << " " << wghts[i]
					<< std::endl;
		}
	}
};
void readdistrs(std::istream &ifs, std::map<std::string, RefDistr> & distrs) {
	std::string line;
	bool indist = false;
	std::string atype;
	double ratio;
	while (ifs.good()) {
		std::getline(ifs, line);
		if (!ifs.good())
			break;
		if (line[0] == '&') {
			std::getline(ifs, line);
			std::string st = line.substr(1);
			std::stringstream sstr(st);
			sstr >> atype;
			sstr >> ratio;
			if (distrs.find(atype) == distrs.end()) {
				distrs.insert(std::make_pair(atype, RefDistr()));
			}
			distrs[atype].rate += ratio;
			distrs[atype].readtimes += 1.0;
			continue;
		}
		std::stringstream sstr(line);
		double r, w;
		sstr >> r >> w;
		int bin = DistBins::distbins().binid(r);
		distrs[atype].wghts[bin] += w;
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

double calcscore(double total_expect, double expect, double observe) {
	double w = 10 / total_expect;
	w = w * w * w;
	double wmin;
	if (observe >= 1 && observe / expect >= 5) {
		observe = observe - 0.8;
	}
	return (observe + w * expect) / ((1.0 + w) * expect);
}
double getscore(double r, std::vector<double> &et) {
	if (r < DistBins::distbins().rmin() || r > DistBins::distbins().rmax())
		return 0.0;
	int idx = DistBins::distbins().binid(r);
	double rc = DistBins::distbins().bincenter(idx);
	int nidx = idx + 1;
	if (r < rc)
		nidx = idx - 1;
//	if (nidx < 0 || nidx >= DistBins::distbins().num_bins())
	return et.at(idx);
//	double w1 = (r - rc) / (DistBins::distbins().bincenter(nidx) - rc);
//	return w1 * (et.at(nidx)) + (1 - w1) * (et.at(idx));
}
double scorelm(const ContactDetails &dtl,
		std::map<std::string, std::vector<double>> &ene) {
	assert(dtl.indirect != 0);
	std::string latype = getlatype(dtl.latype);
	std::string matype = dtl.mligand;
	std::string key = latype + ":" + matype;
	if (ene.find(key) == ene.end())
		return 0.0;
	return getscore(dtl.mdistance, ene[key]);
}
double scorelp(const ContactDetails &dtl,
		std::map<std::string, std::vector<double>> &ene) {
	std::string latype = getlatype(dtl.latype);
	if (dtl.indirect != 0)
		latype = dtl.mligand;
	std::string patype = getpatype(dtl.prname, dtl.paname);
	std::string key = latype + ":" + patype;
	if (ene.find(key) == ene.end())
		return 0.0;
	double r = dtl.distance;
	return getscore(r, ene[key]);
}
void runstatistics(int, char **, std::map<std::string, std::vector<double>> &,
		double &, double);
int main(int argc, char** argv) {
	std::map<std::string, std::vector<double>> ene;
	std::map<std::string, std::vector<double>> ene_new;
	std::map<std::string, std::vector<double>> ene_trial;
	std::map<std::string, std::vector<double>> dene;
	std::ifstream ifsene("ENETable_start.dat");
	while (ifsene.good()) {
		std::string line;
		std::getline(ifsene, line);
		if (!ifsene.good())
			break;
		std::stringstream sstr(line);
		std::string key;
		sstr >> key;
		ene[key] = std::vector<double>(DistBins::distbins().num_bins(), 0.0);
		for (auto &e : ene[key])
			sstr >> e;
	}

	double sumold = 0.0;
	double rate = 0.04;
	double sum = 0.0;
	ene_new = ene;
	runstatistics(argc, argv, ene_new, sumold, rate);
	for (auto &e : ene) {
		dene[e.first] = std::vector<double>();
		for (int n = 0; n < e.second.size(); ++n) {
			dene[e.first].push_back((ene_new[e.first][n] - e.second[n]));
		}
	}
	ene_trial = ene_new;
	double rt = 1.0;
	for (int i = 0; i < 2000; ++i) {
		ene_new = ene_trial;
		runstatistics(argc, argv, ene_new, sum, rate);
		if (sum < sumold) {
			std::ofstream ofsc("ENEChange.dat");
			for (auto &e : ene) {
				ofsc << e.first;
				for (int n = 0; n < e.second.size(); ++n)
					ofsc << " " << ene_trial[e.first][n] - e.second[n];
				ofsc << std::endl;
			}
			ene = ene_trial;
			ene_trial = ene_new;
			if (rate <= 0.01) {
				for (auto & e : ene_trial) {
					if (NSPdstl::RandomEngine<>::getinstance().realrng(0.0, 1.0)()
							> 0.2)
						e.second = ene[e.first];
				}
			}
			sumold = sum;
			rt = 1.0;
			std::ofstream ofs("ENETable.dat");
			for (auto &e : ene) {
				ofs << e.first;
				dene[e.first] = std::vector<double>();
				for (int n = 0; n < e.second.size(); ++n) {
					dene[e.first].push_back(
							(ene_new[e.first][n] - e.second[n]));
				}
				for (auto d : e.second)
					ofs << " " << d;
				ofs << std::endl;
			}
		} else {
//			ene_trial = ene;

			if (rate > 0.01) {
				rt = 0.5 * rt;
				rate = rate * 0.5;
			}
			for (auto & e : ene_trial) {
				if (rate <= 0.01) {
					rt = 1.0;
					if (NSPdstl::RandomEngine<>::getinstance().realrng(0.0, 1.0)()
							> 0.2)
						rt = 0.0;
				}
//				std::cout << rt << std::endl;
				for (int n = 0; n < e.second.size(); ++n) {
					e.second[n] = ene[e.first][n] + rt * dene[e.first][n];
				}
			}
		}
	}
}
double ligandweight(const std::string & molname) {
	static std::map<std::string, double> wghts;
	static bool first { true };
	if (first) {
		first = false;
		std::ifstream ifs("titles.txt");
		std::string line;
		while (ifs.good()) {
			std::getline(ifs, line);
			if (!ifs.good())
				break;
			std::string lname = line.substr(5, 3);
			if (wghts.find(lname) == wghts.end()) {
				wghts[lname] = 0.0;
			}
			wghts[lname] += 1.0;
		}
		for (auto &w : wghts)
			w.second = 1.0 / w.second;
	}
	std::string lname = molname.substr(5, 3);
	if (wghts.find(lname) == wghts.end())
		return -1.0;
	return wghts[lname];
}
void runstatistics(int argc, char **argv,
		std::map<std::string, std::vector<double>> &ene, double &sum2,
		double rate) {
	std::map<std::string, double> num_ligandatoms;
	std::map<std::string, double> num_platoms;
	std::map<std::string, double> num_platoms_w;
	std::map<std::string, std::map<std::string, double>> num_pmatoms;
	std::map<std::string, double> num_mlatoms;
	std::map<std::string, double> num_mlatoms_w;
	std::map<std::string, std::vector<double>> num_lpcontacts;
	std::map<std::string, std::vector<double>> num_lpcontacts_w;
	std::map<std::string, std::vector<double>> num_lmcontacts;
	std::map<std::string, std::vector<double>> num_lmcontacts_w;
	std::map<std::string, std::map<std::string, std::vector<double>>>num_pmcontacts;
	const DistBins &bins = DistBins::distbins();
	std::ifstream ifs((std::string(argv[1])));
	std::map<std::string, AtomType> map = AtomType::getmap();
	std::map<std::string, std::string> nmmap;
	int mode1 = std::atoi(argv[3]);
	int mode2 = std::atoi(argv[4]);
	if (mode1 != 2)
		getlatype = subsitedesign::getcodename0;
	else
		getlatype = subsitedesign::getcodename;
	if (mode2 != 2)
		getpatype = subsitedesign::proteinatomtype;
	else
		getpatype = subsitedesign::proteinatomtype2;
	while (ifs.good()) {
		int nlatoms;
		std::vector<std::string> atypes;
		std::set<int> excluded;
		std::string molname;
		std::vector<ContactDetails> dts = getdtlsnextmol(ifs, molname, nlatoms,
				atypes, excluded);
		if (dts.empty() || nlatoms == 0)
			continue;
		int idx = -1;
		double wlgd = ligandweight(molname);
		if (wlgd < 0)
			break;
		for (auto &at : atypes) {
			++idx;
			std::string latype = getlatype(at);
			nmmap[latype] = subsitedesign::getcodename(at);
			if (latype[0] == 'H')
				continue;
			if (excluded.find(idx) != excluded.end())
				continue;
			if (num_ligandatoms.find(latype) == num_ligandatoms.end())
				num_ligandatoms[latype] = 0.0;
			num_ligandatoms[latype] += wlgd;
		}
		MolContactSum molsum(dts);
		typedef std::pair<int, int> ResidueID;
		std::map<ResidueID, double> prscores;
		for (auto &e : molsum.dcsetsum) {
			if (prscores.find(e.first) == prscores.end())
				prscores[e.first] = 0.0;
			for (auto &sd : e.second) {
				for (auto &c : sd.second)
					prscores[e.first] += scorelp(*c, ene);
			}
		}
		for (auto &prs : prscores) {
//			std::cout <<prs.second<<std::endl;
			prs.second = exp(prs.second);
		}
		for (auto &e : molsum.dcsetsum) {
			std::string prname = e.second.begin()->second[0]->prname;
			double wr = prscores[e.first];
			for (auto &sd : e.second) {
				std::string patname = getpatype(prname, sd.first);
				if (patname.empty())
					continue;
				if (num_platoms.find(patname) == num_platoms.end()) {
					num_platoms[patname] = 0.0;
					num_platoms_w[patname] = 0.0;
				}
				num_platoms[patname] += wlgd;
				num_platoms_w[patname] += wr * wlgd;
				for (auto d : sd.second) {
					std::string latype;
					try {
						latype = getlatype(d->latype);
					} catch (std::exception &e) {
						std::cout << d->latype << std::endl;
						exit(1);
					}
					std::string con_type = latype + ":" + patname;
					if (num_lpcontacts.find(con_type) == num_lpcontacts.end()) {
						num_lpcontacts[con_type] = std::vector<double>(
								bins.num_bins(), 0.0);
						num_lpcontacts_w[con_type] = std::vector<double>(
								bins.num_bins(), 0.0);
					}
					int idx = bins.binid(d->distance);
					num_lpcontacts[con_type][idx] += wlgd;
					num_lpcontacts_w[con_type][idx] += wr * wlgd;
				}
			}
		} //direct contacts
		std::map<ResidueID, double> mrscores;
		for (auto &me : molsum.indcsetsums) {
			if (mrscores.find(me.first) == mrscores.end())
				mrscores[me.first] = 0.0;
			mrscores[me.first] += scorelm(
					*(me.second.begin()->second.begin()->second[0]), ene);
		}
		for (auto &mrs : mrscores)
			mrs.second = exp(mrs.second);
		for (auto &me : molsum.indcsetsums) {
			std::string mname =
					me.second.begin()->second.begin()->second[0]->mligand;
			if (num_mlatoms.find(mname) == num_mlatoms.end()) {
				num_mlatoms[mname] = 0.0;
				num_mlatoms_w[mname] = 0.0;
			}
			double wm = mrscores[me.first];
			num_mlatoms[mname] += wlgd;
			num_mlatoms_w[mname] += wm * wlgd;
			for (auto &id : molsum.inddistances[me.first]) {
				/*				if(id.first >=atypes.size()){
				 break;
				 }
				 if(map.find(atypes[id.first])==map.end()) {
				 break;
				 }*/
				std::string latype;
				try {
					latype = getlatype(atypes[id.first]);
				} catch (std::exception &e) {
					std::cout << atypes.size() << " " << id.first << std::endl;
					for (auto &t : atypes)
						std::cout << " " << t;
					std::cout << std::endl;
					exit(1);
				}
				std::string con_type = latype + ":" + mname;
				if (num_lmcontacts.find(con_type) == num_lmcontacts.end()) {
					num_lmcontacts[con_type] = std::vector<double>(
							bins.num_bins(), 0.0);
					num_lmcontacts_w[con_type] = std::vector<double>(
							bins.num_bins(), 0.0);
				}
				num_lmcontacts[con_type][bins.binid(id.second)] += wlgd;
				num_lmcontacts_w[con_type][bins.binid(id.second)] += wm * wlgd;
			}
			if (num_pmatoms.find(mname) == num_pmatoms.end()) {
				num_pmatoms[mname] = std::map<std::string, double>();
			}
			if (num_pmcontacts.find(mname) == num_pmcontacts.end()) {
				num_pmcontacts[mname] = std::map<std::string,
						std::vector<double>>();
			}
			for (auto &e : me.second) {
				std::string prname = e.second.begin()->second[0]->prname;
				for (auto &sd : e.second) {
					std::string patname = getpatype(prname, sd.first);
					if (patname.empty())
						continue;
					if (num_pmatoms[mname].find(patname)
							== num_pmatoms[mname].end()) {
						num_pmatoms[mname][patname] = 0.0;
					}
					num_pmatoms[mname][patname] += wlgd;
					for (auto d : sd.second) {
						int idx = bins.binid(d->distance);
						std::string ct = patname;
						if (num_pmcontacts[mname].find(ct)
								== num_pmcontacts[mname].end()) {
							num_pmcontacts[mname][ct] = std::vector<double>(
									bins.num_bins(), 0.0);
						}
						num_pmcontacts[mname][ct][idx] += wlgd;
					}
				}
			}
		} //indirect
	} //ifs good
	double ntot_latoms = 0;
	for (auto &n : num_ligandatoms)
		ntot_latoms += n.second;
	double ntot_platoms = 0;
	for (auto &n : num_platoms)
		ntot_platoms += n.second;
	double ntot_mlatoms = 0;
	for (auto &n : num_mlatoms)
		ntot_mlatoms += n.second;
	double ntot_platoms_w = 0;
	for (auto &n : num_platoms_w)
		ntot_platoms_w += n.second;
	double ntot_mlatoms_w = 0;
	for (auto &n : num_mlatoms_w)
		ntot_mlatoms_w += n.second;

	double ntot_lpcontacts = 0;
	for (auto &e : num_lpcontacts) {
		for (auto &d : e.second)
			ntot_lpcontacts += d;
	}
	double ntot_lmcontacts = 0;
	for (auto &e : num_lmcontacts) {
		for (auto &d : e.second)
			ntot_lmcontacts += d;
	}
	double ntot_lpcontacts_w = 0;
	for (auto &e : num_lpcontacts_w) {
		for (auto &d : e.second)
			ntot_lpcontacts_w += d;
	}
	double ntot_lmcontacts_w = 0;
	for (auto &e : num_lmcontacts_w) {
		for (auto &d : e.second)
			ntot_lmcontacts_w += d;
	}

	double ntot_lcontacts = ntot_lpcontacts + ntot_lmcontacts;
	double ntot_pmlatoms = ntot_platoms + ntot_mlatoms;
	double ntot_lcontacts_w = ntot_lpcontacts_w + ntot_lmcontacts_w;
	double ntot_pmlatoms_w = ntot_platoms_w + ntot_mlatoms_w;
	double rate_av=ntot_lcontacts/ntot_latoms;
	double rate_av_w=ntot_lcontacts_w/ntot_latoms;
	double rate_avl = ntot_lpcontacts / ntot_latoms;
	double rate_avl_w = ntot_lpcontacts_w / ntot_latoms;
	double rate_avm = ntot_lmcontacts/ntot_latoms;
	double rate_avm_w=ntot_lmcontacts_w/ntot_latoms;
	double nmetal=ntot_mlatoms-num_mlatoms["HOH"];
//	std::map<std::string,double> refratio;
//	std::map<std::string,std::vector<double>> refwghts;
	//for test
	/*	for(auto & e:num_ligandatoms){
	 refratio[e.first]=1.0;
	 refwghts[e.first]=std::vector<double>(bins.num_bins(),1.0);
	 }*/

	std::ifstream ifsd((std::string(argv[2])));
	std::map<std::string, RefDistr> distrs;
	readdistrs(ifsd, distrs);
	for (auto &d : distrs) {
		d.second.average();
		/*		std::cout <<"&"<<std::endl;
		 std::cout <<"#"<<d.first;
		 d.second.print(std::cout);*/
	}
//	exit(0);
	std::string mtype_equiv = "O.31";
	sum2 = 0.0;
	std::string outfile = std::string("ScoreTable_") + std::string(argv[3])
			+ "_" + std::string(argv[4]) + ".dat";
	std::ofstream ofscst(outfile);
	for (auto & e : num_lpcontacts_w) {
		std::string latype = e.first.substr(0, e.first.find(':'));
		std::string patype = e.first.substr(e.first.find(':') + 1);
		if (ene.find(e.first) == ene.end()) {
			ene[e.first] = std::vector<double>(DistBins::distbins().num_bins(),
					0.0);
		}
		std::vector<double> &et = ene[e.first];
		double nc = 0;
		double nc_w = 0.0;
		for (auto n : num_lpcontacts.at(e.first))
			nc += n;
		for (auto n : e.second) {
			nc_w += n;
		}
		double nc_expect_w = rate_av_w * num_ligandatoms[latype]
				* num_platoms_w[patype] / ntot_pmlatoms_w
				* distrs[nmmap[latype]].rate;
		double nc_expect = rate_av * num_ligandatoms[latype]
				* num_platoms[patype] / ntot_pmlatoms
				* distrs[nmmap[latype]].rate;
//		std::vector<double> score(bins.num_bins(),0);
		ofscst << "&" << std::endl;
		ofscst << latype << " " << patype << " " << nc << " " << nc_expect
				<< std::endl;
		for (int i = 0; i < bins.num_bins(); ++i) {
			double expect = nc_expect * distrs[nmmap[latype]].wghts[i] + 1.e-5;
			ofscst << bins.bincenter(i) << " " << et[i] << " "
					<< num_lpcontacts[e.first][i] << " " << expect << std::endl;
			double expect_w = nc_expect_w * distrs[nmmap[latype]].wghts[i]
					+ 1.e-5;
			double nscale=nc>nc_expect? nc:nc_expect;
			double s = calcscore(nscale, expect_w, e.second[i]);
			s = log(s);
			double rt = rate;
			if (e.second[i] != 0.0) {
				sum2 += s * s;
			}
			double eb=-log(calcscore(nscale,expect,num_lpcontacts[e.first][i]));
			double emin,emax;
			if(eb>0.0){
				emin=0.0;
				emax=eb;
			} else {
				emin=eb;
				emax=0.0;
			}
			et[i] = (1 - rt) * et[i] - rt * s;
			if(et[i]<emin) et[i]=emin;
			if(et[i]>emax) et[i]=emax;
		}
	}
	for (auto & e : num_lmcontacts_w) {
		std::string latype = e.first.substr(0, e.first.find(':'));
		std::string matype = e.first.substr(e.first.find(':') + 1);
		if (ene.find(e.first) == ene.end())
			ene[e.first] = std::vector<double>(DistBins::distbins().num_bins(),
					0.0);
		std::vector<double> &et = ene[e.first];
		double nc_w = 0;
		double nc = 0;
		for (auto n : num_lmcontacts.at(e.first))
			nc += n;
		for (auto n : e.second) {
			nc_w += n;
		}
		/*		double nc_expect = rate_av * num_ligandatoms[latype]
		 * num_mlatoms[matype] / ntot_pmlatoms
		 * distrs[nmmap[latype]].rate;
		 double nc_expect_w = rate_av_w * num_ligandatoms[latype]
		 * num_mlatoms_w[matype] / ntot_pmlatoms_w
		 * distrs[nmmap[latype]].rate;*/
		double nc_expect = rate_av * num_ligandatoms[latype]*
							ntot_mlatoms/ntot_pmlatoms* distrs[nmmap[latype]].rate;
		double nc_expect_w = rate_av_w * num_ligandatoms[latype]*
				ntot_mlatoms/ntot_pmlatoms
				* distrs[nmmap[latype]].rate;
//		if(matype !="HOH"){
			nc_expect=nc_expect*num_mlatoms[matype]/ntot_mlatoms;
			nc_expect_w=nc_expect_w*num_mlatoms[matype]/ntot_mlatoms;
//		}
//		std::vector<double> score(bins.num_bins(),0);

		ofscst << "&" << std::endl;
		ofscst << latype << " " << matype << " " << nc << " " << nc_expect
				<< std::endl;
		for (int i = 0; i < bins.num_bins(); ++i) {
			double expect = nc_expect * distrs[nmmap[latype]].wghts[i] + 1.e-5;
			ofscst << bins.bincenter(i) << " " << et[i] << " "
					<< num_lmcontacts[e.first][i] << " " << expect << std::endl;
			double expect_w = nc_expect_w * distrs[nmmap[latype]].wghts[i]
					+ 1.e-5;
			double nscale=nc>nc_expect?nc:nc_expect;
			double s = calcscore(nscale, expect_w, e.second[i]);
			s = log(s);
			double rt = rate;
			if (e.second[i] != 0) {
				sum2 += s * s;
			}
			et[i] = (1 - rt) * et[i] - rt * s;
			double eb=-log(calcscore(nscale,expect,num_lmcontacts[e.first][i]));
			double emin,emax;
			if(eb>0.0){
				emin=0.0;
				emax=eb;
			} else {
				emin=eb;
				emax=0.0;
			}
			et[i] = (1 - rt) * et[i] - rt * s;
			if(et[i]<emin) et[i]=emin;
			if(et[i]>emax) et[i]=emax;

//			std::cout <<bins.bincenter(i) <<" "<<score[i]
//					<<" "<<e.second[i]<<" "<<expect<<std::endl;
		}
	}
	std::cout << "Sum of ediff2: " << sum2 << std::endl;
	std::map<std::string, double> ref_pref;
	for (auto & npl : num_platoms) {
		ref_pref[npl.first] = npl.second / ntot_platoms; //protein atomtype distr, ligand
	}
	for (auto &me : num_pmcontacts) {
		std::string matype = me.first;
		double ntot_pm = 0.0;
		for (auto &n : me.second) {
			for (auto d : n.second)
				ntot_pm += d;
		}
		double ntot_pmatoms = 0.0;
		for (auto &n : num_pmatoms[me.first]) {
			ntot_pmatoms += n.second;
		}
		for (auto &pe : me.second) {
			std::string patype = pe.first;
			std::vector<double> score(bins.num_bins(), 0.0);
			double pref = num_pmatoms[matype][patype] / ntot_pmatoms;
			double ncon = 0.0;
			for (auto d : pe.second)
				ncon += d;
			ofscst << "&" << std::endl;
			ofscst << matype << " " << patype << " " << ncon << " " << ncon
					<< std::endl;
			for (int i = 0; i < bins.num_bins(); ++i) {
				double expect = ncon * distrs[mtype_equiv].wghts[i] + 1.e-5;
				score[i] = pref / ref_pref[patype] * (pe.second[i] + 1.e-5)
						/ expect;
//						* calcscore(ncon, expect, pe.second[i]);
				ofscst << bins.bincenter(i) << " " << -log(score[i]) << " "
						<< pe.second[i] << " " << expect << std::endl;
			}
		}
	}
}
