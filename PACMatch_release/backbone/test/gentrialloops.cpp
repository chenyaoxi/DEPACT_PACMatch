/*
 * gentrialloops.cpp
 *
 *  Created on: 2017年4月29日
 *      Author: hyliu
 */

#include "backbone/mainchain.h"
#include "backbone/energyfunctions.h"
#include "backbone/closingloop.h"
#include "dstl/randomengine.h"
#include "dstl/permutation.h"
#include "dstl/topn.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <memory>
using namespace NSPproteinrep;

bool readsegments(std::istream &is,
		std::vector<std::vector<BackBoneSite>> *segments) {
	int nsegments;
	std::istringstream iss;
	std::vector<int> segmentsizes;
	char buffer[120];
	is.getline(buffer, 120);
	iss.str(std::string(buffer));
	iss >> nsegments;
	segmentsizes.resize(nsegments, 0);
	std::vector<std::pair<int,int>> enddel(nsegments,{0,0});
	for (int i = 0; i < nsegments; ++i) {
		is.getline(buffer, 120);
		iss.str(std::string(buffer));
		iss.seekg(0);
		iss >> segmentsizes[i] >>enddel[i].first >>enddel[i].second;
	}
	bool success = false;
	for (int i = 0; i < nsegments; ++i) {
		segments->push_back(std::vector<BackBoneSite>());
		std::vector<BackBoneSite> tmp;
		success = readbackbonesites(is, segmentsizes[i], tmp);
		int newsize=segmentsizes[i]-enddel[i].first-enddel[i].second;
		segments->back().resize(newsize);
		std::copy(tmp.begin()+enddel[i].first,
				tmp.begin()+enddel[i].first+newsize,segments->back().begin());
		if (!success)
			return success;

	}
	return success;
}

int loops1s2l(const std::vector<std::vector<BackBoneSite>> & segments, int s1,
		int s2, int length, int nstop,
		std::vector<std::shared_ptr<std::vector<BackBoneSite>>>*solutions,
		std::vector<double> *enes) {
	std::vector<int> order, looplengths;
	order.push_back(s1);
	order.push_back(s2);
	for (int i = 0; i < segments.size(); ++i) {
		if (s1 != i && s2 != i)
		order.push_back(i);
	}
	looplengths.push_back(length);
	for (int i = 1; i < segments.size() - 1; ++i) {
		looplengths.push_back(7);
	}
	MainChain mc;
	mainchainfromsegments(segments, order, looplengths, &mc);
	int loopstart = segments[order[0]].size() - 1;
	int loopend = loopstart + looplengths[0] + 2;
	ClosingLoop cl(mc, loopstart, loopend, ClosingLoop::NEWSEQ_NEWCONF);
	EnergyTerms eterms;
	eterms.addstericclash(10.0);
	eterms.addtorsionene();
	EnergyTerms eterms_rank;
	eterms_rank.addstericclash(10.0);
//	eterms_rank.addtorsionvecene();
	eterms_rank.addtetrasefene(0.17);
	int ntotalsol=0;
	NSPdstl::TopN<std::shared_ptr<std::vector<BackBoneSite>>> topn(nstop);
	for (int ntry = 0; ntry < 1000*nstop; ++ntry) {
		int nsol = cl.solve();
		for (int i = 0; i < nsol; ++i) {
//			if(ntotalsol >=nstop) break;
			std::shared_ptr<std::vector<BackBoneSite>> newloopi(
			new std::vector<BackBoneSite>());
			cl.getsitessolution(i, newloopi.get());
			std::map<std::string, double> energies;
			double eclash = calcloopenergy(mc, loopstart, newloopi->size(),
			*newloopi, eterms, energies, eterms.emaxallowed);
//			if (eclash >= eterms.emaxallowed)
//				std::cout <<eclash <<std::endl;
			if(eclash >= eterms.emaxallowed || energies["TorsionEne"] >-1.0*(double) (newloopi->size())) continue;
			double erank=calcloopenergy(mc, loopstart, newloopi->size(),
					*newloopi, eterms_rank, energies, eterms_rank.emaxallowed);
			topn.push(newloopi,erank);
			++ntotalsol;
		}
		if(ntry>=5*nstop && ntotalsol >=nstop) break;
	}
	std::vector<std::pair<std::shared_ptr<std::vector<BackBoneSite>>,double>>
	      res=topN2vector(topn);
	std::cout <<"\tLoop energies:"<< std::endl;
	for (auto &r:res) {
		solutions->push_back(r.first);
		enes->push_back(r.second);
		std::cout <<"\t" << r.second;
	}
	std::cout <<std::endl;
	return res.size();
}
int loops1s2(const std::vector<std::vector<BackBoneSite>> & segments, int s1,
		int s2, int minlength, int maxlength, int nstop,
		std::vector<std::shared_ptr<std::vector<BackBoneSite>>>*solutions,
		std::vector<double> *enes) {
	int nsol=0;
	for(int len=minlength; len<=maxlength;++len) {
		int nsol_len =loops1s2l(segments,s1,s2,len,nstop,solutions,enes);
		std::cout << "\tNumber of loops of length " << len <<": " <<nsol_len <<std::endl;
		nsol +=nsol_len;
		if(len ==maxlength && nsol < nstop && maxlength<10) ++maxlength;
//		if(nsol >= nstop) break;
	}
	return nsol;
}
typedef std::shared_ptr<std::vector<BackBoneSite>> LoopPtr;
typedef std::vector<LoopPtr> OneLoopSolutions;
void buildtrialloops(const std::vector<std::vector<BackBoneSite>> &segments,
		std::map<std::pair<int, int>, OneLoopSolutions> *solutions,
		std::map<std::pair<int, int>, std::vector<double>> *enes) {
	int minlength = 1;
	int maxlength = 5;
	int nstop = 15;
	for (int s1 = 0; s1 < segments.size(); ++s1) {
		for (int s2 = 0; s2 < segments.size(); ++s2) {
			if (s2 == s1)
				continue;
			std::cout << "Trying to build " << s1 << "-" << s2 << " loops... "<<std::endl;
			std::pair<int, int> index = std::make_pair(s1, s2);
			solutions->insert(std::make_pair(index, OneLoopSolutions()));
			enes->insert(std::make_pair(index, std::vector<double>()));
			int nsol = loops1s2(segments, s1, s2, minlength, maxlength, nstop,
					&(solutions->at(index)), &(enes->at(index)));
			std::cout << nsol << " clash-free loops built." << std::endl;
		}
	}
}
void tryends(const std::vector<std::vector<BackBoneSite>> &segments) {
	int minlength = 1;
	int maxlength = 5;
	int nstop = 15;
	std::vector<std::pair<int,int>> enddels{{0,0},{1,0},{2,0},{3,0},
		{0,1},{1,1},{2,1},{3,1},{0,2},{1,2},{2,2},{3,2},{0,3},
		{1,3},{2,3},{3,3}
	};
	for (int s1 = 0; s1 < segments.size(); ++s1) {
		for (int s2 = 0; s2 < segments.size(); ++s2) {
			if (s2 == s1)
				continue;
			std::cout << "Trying to build " << s1 << "-" << s2 << " loops... "<<std::endl;
			for(auto &ed:enddels) {
				std::vector<std::vector<BackBoneSite>> newsgmnts;
				newsgmnts.resize(segments.size());
				bool dotry=true;
				for(int m=0; m<segments.size(); ++m) {
					if(m==s1) {
						int isize=segments[m].size()-ed.first;
						if(isize <=1) {dotry=false; break;}
						newsgmnts[m].resize(isize);
						std::copy(segments[m].begin(),segments[m].begin()+isize,
								newsgmnts[m].begin());
					} else if(m==s2){
						int isize=segments[m].size()-ed.second;
						if(isize<=1) {dotry=false;break;}
						newsgmnts[m].resize(isize);
							std::copy(segments[m].begin()+ed.second,
									segments[m].end(),
									newsgmnts[m].begin());
					} else {
						newsgmnts[m]=segments[m];
					}
				}
				if(!dotry) continue;
				std::pair<int, int> index = std::make_pair(s1, s2);
				OneLoopSolutions sl;
				std::vector<double> ene;
				int nsol = loops1s2(newsgmnts, s1, s2, minlength, maxlength, nstop,
					&sl, &ene);
				if(nsol>0) {
					std::cout <<"For S"<<s1<<"-"<<"S"<<s2<<"loops, ";
					std::cout <<"solutions exists for enddels: "<< ed.first <<"\t"<<ed.second<<std::endl;
				}
			}
		}
	}
}
void selectnextloop(const std::vector<int> & segmentorder,
		const std::map<std::pair<int, int>, OneLoopSolutions> & candidateloops,
		const  std::map<std::pair<int, int>, std::vector<double>> &candidateenergies,
		std::vector<int> & selectedloops,
		std::vector<std::pair<std::vector<LoopPtr>,double>> *result) {
	int nselected = selectedloops.size();
	if (nselected == segmentorder.size() - 1) {
		result->push_back(std::make_pair(std::vector<LoopPtr>(),0.0));
		std::vector<LoopPtr> &res = result->back().first;
		double & etotal=result->back().second;
		for (int sl = 0; sl < nselected; sl++) {
			int ss1 = segmentorder[sl];
			int ss2 = segmentorder[sl + 1];
			res.push_back(
					candidateloops.at(std::make_pair(ss1, ss2))[selectedloops[sl]]);
			etotal += candidateenergies.at(std::make_pair(ss1, ss2))[selectedloops[sl]];
		}
		return;
	}
	int s1 = segmentorder[nselected];
	int s2 = segmentorder[nselected + 1];
	const OneLoopSolutions &candidates = candidateloops.at(
			std::make_pair(s1, s2));
	StericClash *efunc = &(eneterminstance<StericClash>());
	int max_used_can=1;
	std::map<int,int> nused_can;
	for (int ican = 0; ican < candidates.size(); ++ican) {
		bool accept_ican = true;
		const std::vector<BackBoneSite> *loopcan = candidates[ican].get();
		int len=loopcan->size();
		if(nused_can.find(len) != nused_can.end()) {
			if(nused_can[len]>=max_used_can) continue;
		}
		for (int sl = 0; sl < nselected; sl++) {
			int ss1 = segmentorder[sl];
			int ss2 = segmentorder[sl + 1];
			const std::vector<BackBoneSite> *loopsl = candidateloops.at(
					std::make_pair(ss1, ss2))[selectedloops[sl]].get();
			double e = looploopenergy(*loopcan, *loopsl, efunc,
					efunc->eclash());
			if (e >= efunc->eclash()) {
				accept_ican = false;
			}
		}
		if (accept_ican) {
			selectedloops.push_back(ican);
			int oldnres=result->size();
			selectnextloop(segmentorder, candidateloops, candidateenergies,selectedloops, result);
			int nres=result->size() - oldnres;
			if(nres >0) {
				if(nused_can.find(len) == nused_can.end()){
					nused_can.insert(std::make_pair(len,0));
				}
				nused_can[len] +=(result->size() - oldnres);
			}
			selectedloops.pop_back();
		}
	}
}

int selectloops(const std::vector<int> & segmentorder,
		const std::map<std::pair<int, int>, OneLoopSolutions> & candidateloops,
		const  std::map<std::pair<int, int>, std::vector<double>> &candidateenergies,
		std::vector<std::pair<std::vector<LoopPtr>,double>> *result) {
	bool possible = true;
	for (int i = 0; i < segmentorder.size() - 1; ++i) {
		if (candidateloops.at(
				std::make_pair(segmentorder[i], segmentorder[i + 1])).empty()) {
			possible = false;
			return 0;
		}
	}
	if (possible) {
		std::vector<int> selectedloops;
		selectnextloop(segmentorder, candidateloops, candidateenergies,selectedloops, result);
	}
	return result->size();
}

std::string namelinkedchain(const std::vector<int> & segmentorder,
		const std::vector<LoopPtr> & linkers){
		std::string str;
		for (int i=0; i<segmentorder.size()-1;++i) {
			str =str + "S" + std::to_string(segmentorder[i]);
			str =str + "_L"+ std::to_string(linkers[i]->size()-2) +"_";
		}
		str = str+"S"+std::to_string(segmentorder.back());
		return str;
}

void buildlinkedchain(const std::vector<std::vector<BackBoneSite>> & segments,
		const std::vector<int> & segmentorder,
				const std::pair<std::vector<LoopPtr>,double> & linkers,
				MainChain *mc,double *eneperloopresidue){
		std::vector<int> looplens;
		int nloopresidues=0;
		for (auto &lp:linkers.first) {
			looplens.push_back(lp->size()-2);
			nloopresidues += looplens.back()+2;
		}
		*eneperloopresidue=linkers.second/(double) nloopresidues;
		mainchainfromsegments(segments, segmentorder, looplens, mc);
		int posi=0;
		for(int i=0;i<linkers.first.size();++i) {
			posi +=segments[segmentorder[i]].size()-1;
			std::copy(linkers.first[i]->begin(),linkers.first[i]->end(),
					mc->begin()+posi);
			posi += linkers.first[i]->size()-1;
		}
		mc->resetresseq();
}
int main(int argc, char **argv) {
	std::string filename(argv[1]);
	std::ifstream ifs;
	ifs.open(filename.c_str());
	std::vector<std::vector<BackBoneSite>> segments;
	readsegments(ifs, &segments);
	if(argc>=3){
		if(std::string(argv[2]) == "trydelends"){
			tryends(segments);
			exit(0);
		}
	}
	typedef std::vector<std::shared_ptr<std::vector<BackBoneSite>>>OneLoopSolutions;
	std::map<std::pair<int, int>, OneLoopSolutions> candidateloops;
	std::map<std::pair<int, int>,std::vector<double>>candidateenergies;
	buildtrialloops(segments, &candidateloops,&candidateenergies);
	std::vector<int> origin_order;
	for (unsigned int i = 0; i < segments.size(); ++i) {
		origin_order.push_back(i);
	}
	int np = NSPdstl::Permutation<int>::factorial(origin_order.size());

	for (unsigned int k = 0; k < np; ++k) {
		std::vector<int> order = NSPdstl::Permutation<int>::getPermutation(
				origin_order, k);
		std::cout << "permutation " << k << ": ";
		for (auto o : order)
			std::cout << "\t" << o;
		std::cout << std::endl;
		std::vector<std::pair<std::vector<LoopPtr>,double>> result;
		int ncomplete = selectloops(order, candidateloops, candidateenergies, &result);
		std::cout << " Number of completely-linked solutions: " << ncomplete
				<< std::endl;
		for(int i=0;i<ncomplete; ++i) {
			std::string chainname=namelinkedchain(order,result[i].first);
			MainChain mc;
			double eneperloopresidue;
			buildlinkedchain(segments,
					order,
					result[i], &mc,&eneperloopresidue);
		     std::cout <<"Chain " <<chainname <<" energy per loop resdiue: "<<
		    		 eneperloopresidue<<std::endl;
			 std::ofstream ofs;
			 std::string filename=chainname+".pdb";
			 ofs.open(filename.c_str());
			 writeSitesToPDB(ofs, mc);
			 ofs.close();
			 filename=chainname+".dat";
			 ofs.open(filename.c_str());
			 mc.write(ofs);
		     ofs.close();
		}
	}
	/*	for (int s1 = 0; s1 < segments.size(); ++s1) {
	 for (int s2 = 0; s2 < segments.size(); ++s2) {
	 if (s2 == s1)
	 continue;
	 std::pair<int, int> index = std::make_pair(s1, s2);
	 auto &osl = solutions[index];
	 std::cout << s1 << "-" << s2 << " loops: " << osl.size()
	 << std::endl;
	 for (auto & s : osl) {
	 std::cout << "\tlooplength: " << s->size() - 2 << std::endl;
	 }
	 }
	 }
	 */
	/*	int i = 0;
	 for (auto iter = solutions.begin(); iter != solutions.end(); ++iter) {
	 std::copy((*iter)->begin(), (*iter)->end(), mc.begin() + loopstart);
	 mc.resetresseq();
	 filename = "temp" + std::to_string(i++) + ".pdb";
	 std::ofstream ofs;
	 ofs.open(filename.c_str());
	 writeSitesToPDB(ofs, mc);
	 }
	 */
}

