/*
 * aminoacidseq.cpp
 *
 *  Created on: 2017年6月16日
 *      Author: hyliu
 */

#include "proteinrep/aminoacidseq.h"
#include <cassert>
using namespace NSPproteinrep;

std::map<char,std::string> AminoAcidSeq::CodeName{
	{'A',"ALA"},{'C',"CYS"},{'D',"ASP"},{'E',"GLU"},{'F',"PHE"},{'G',"GLY"},
	{'H',"HIS"},{'I',"ILE"},{'K',"LYS"},{'L',"LEU"},{'M',"MET"},{'N',"ASN"},
	{'P',"PRO"},{'Q',"ASN"},{'R',"ARG"},{'S',"SER"},{'T',"THR"},{'V',"VAL"},
	{'W',"TRP"},{'Y',"TYR"},{'X',"ANY"}};
std::map<std::string,char> AminoAcidSeq::NameCode{
	{"ALA",'A'},{"CYS",'C'},{"ASP",'D'},{"GLU",'E'},{"PHE",'F'},{"GLY",'G'},
	{"HIS",'H'},{"ILE",'I'},{"LYS",'K'},{"LEU",'L'},{"MET",'M'},{"ASN",'N'},
	{"PRO",'P'},{"ASN",'Q'},{"ARG",'R'},{"SER",'S'},{"THR",'T'},{"VAL",'V'},
	{"TRP",'W'},{"TYR",'Y'}};
std::string AminoAcidSeq::name2code(const std::vector<std::string> &nameseq){
	std::string res;
	for(auto &name:nameseq) res.push_back(name2code(name));
	return res;
}
std::vector<std::string> AminoAcidSeq::code2name(const std::string &codeseq){
	std::vector<std::string> res;
	for(auto c:codeseq) res.push_back(code2name(c));
	return res;
}

std::map<std::pair<char,char>,double> SeqAlignment::blosum62matrix{
	{{'A','A'},4}, {{'A','C'},0}, {{'A','D'},-2}, {{'A','E'},-1}, {{'A','F'},-2},
	{{'A','G'},0}, {{'A','H'},-2}, {{'A','I'},-1}, {{'A','K'},-1}, {{'A','L'},-1},
	{{'A','M'},-1}, {{'A','N'},-2}, {{'A','P'},-1}, {{'A','Q'},-1}, {{'A','R'},-1},
	{{'A','S'},1}, {{'A','T'},0}, {{'A','V'},0}, {{'A','W'},-3}, {{'A','Y'},-2},
	{{'C','A'},0}, {{'C','C'},9}, {{'C','D'},-3}, {{'C','E'},-4}, {{'C','F'},-2},
	{{'C','G'},-3}, {{'C','H'},-3}, {{'C','I'},-1}, {{'C','K'},-3}, {{'C','L'},-1},
	{{'C','M'},-1}, {{'C','N'},-3}, {{'C','P'},-3}, {{'C','Q'},-3}, {{'C','R'},-3},
	{{'C','S'},-1}, {{'C','T'},-1}, {{'C','V'},-1}, {{'C','W'},-2}, {{'C','Y'},-2},
	{{'D','A'},-2}, {{'D','C'},-3}, {{'D','D'},6}, {{'D','E'},2}, {{'D','F'},-3},
	{{'D','G'},-1}, {{'D','H'},-1}, {{'D','I'},-3}, {{'D','K'},-1}, {{'D','L'},-4},
	{{'D','M'},-3}, {{'D','N'},1}, {{'D','P'},-1}, {{'D','Q'},0}, {{'D','R'},-2},
	{{'D','S'},0}, {{'D','T'},-1}, {{'D','V'},-3}, {{'D','W'},-4}, {{'D','Y'},-3},
	{{'E','A'},-1}, {{'E','C'},-4}, {{'E','D'},2}, {{'E','E'},5}, {{'E','F'},-3},
	{{'E','G'},-2}, {{'E','H'},0}, {{'E','I'},-3}, {{'E','K'},1}, {{'E','L'},-3},
	{{'E','M'},-2}, {{'E','N'},0}, {{'E','P'},-1}, {{'E','Q'},2}, {{'E','R'},0},
	{{'E','S'},0}, {{'E','T'},-1}, {{'E','V'},-2}, {{'E','W'},-3}, {{'E','Y'},-2},
	{{'F','A'},-2}, {{'F','C'},-2}, {{'F','D'},-3}, {{'F','E'},-3}, {{'F','F'},6},
	{{'F','G'},-3}, {{'F','H'},-1}, {{'F','I'},0}, {{'F','K'},-3}, {{'F','L'},0},
	{{'F','M'},0}, {{'F','N'},-3}, {{'F','P'},-4}, {{'F','Q'},-3}, {{'F','R'},-3},
	{{'F','S'},-2}, {{'F','T'},-2}, {{'F','V'},-1}, {{'F','W'},1}, {{'F','Y'},3},
	{{'G','A'},0}, {{'G','C'},-3}, {{'G','D'},-1}, {{'G','E'},-2}, {{'G','F'},-3},
	{{'G','G'},6}, {{'G','H'},-2}, {{'G','I'},-4}, {{'G','K'},-2}, {{'G','L'},-4},
	{{'G','M'},-3}, {{'G','N'},0}, {{'G','P'},-2}, {{'G','Q'},-2}, {{'G','R'},-2},
	{{'G','S'},0}, {{'G','T'},-2}, {{'G','V'},-3}, {{'G','W'},-2}, {{'G','Y'},-3},
	{{'H','A'},-2}, {{'H','C'},-3}, {{'H','D'},-1}, {{'H','E'},0}, {{'H','F'},-1},
	{{'H','G'},-2}, {{'H','H'},8}, {{'H','I'},-3}, {{'H','K'},-1}, {{'H','L'},-3},
	{{'H','M'},-2}, {{'H','N'},1}, {{'H','P'},-2}, {{'H','Q'},0}, {{'H','R'},0},
	{{'H','S'},-1}, {{'H','T'},-2}, {{'H','V'},-3}, {{'H','W'},-2}, {{'H','Y'},2},
	{{'I','A'},-1}, {{'I','C'},-1}, {{'I','D'},-3}, {{'I','E'},-3}, {{'I','F'},0},
	{{'I','G'},-4}, {{'I','H'},-3}, {{'I','I'},4}, {{'I','K'},-3}, {{'I','L'},2},
	{{'I','M'},1}, {{'I','N'},-3}, {{'I','P'},-3}, {{'I','Q'},-3}, {{'I','R'},-3},
	{{'I','S'},-2}, {{'I','T'},-1}, {{'I','V'},3}, {{'I','W'},-3}, {{'I','Y'},-1},
	{{'K','A'},-1}, {{'K','C'},-3}, {{'K','D'},-1}, {{'K','E'},1}, {{'K','F'},-3},
	{{'K','G'},-2}, {{'K','H'},-1}, {{'K','I'},-3}, {{'K','K'},5}, {{'K','L'},-2},
	{{'K','M'},-1}, {{'K','N'},0}, {{'K','P'},-1}, {{'K','Q'},1}, {{'K','R'},2},
	{{'K','S'},0}, {{'K','T'},-1}, {{'K','V'},-2}, {{'K','W'},-3}, {{'K','Y'},-2},
	{{'L','A'},-1}, {{'L','C'},-1}, {{'L','D'},-4}, {{'L','E'},-3}, {{'L','F'},0},
	{{'L','G'},-4}, {{'L','H'},-3}, {{'L','I'},2}, {{'L','K'},-2}, {{'L','L'},4},
	{{'L','M'},2}, {{'L','N'},-3}, {{'L','P'},-3}, {{'L','Q'},-2}, {{'L','R'},-2},
	{{'L','S'},-2}, {{'L','T'},-1}, {{'L','V'},1}, {{'L','W'},-2}, {{'L','Y'},-1},
	{{'M','A'},-1}, {{'M','C'},-1}, {{'M','D'},-3}, {{'M','E'},-2}, {{'M','F'},0},
	{{'M','G'},-3}, {{'M','H'},-2}, {{'M','I'},1}, {{'M','K'},-1}, {{'M','L'},2},
	{{'M','M'},5}, {{'M','N'},-2}, {{'M','P'},-2}, {{'M','Q'},0}, {{'M','R'},-1},
	{{'M','S'},-1}, {{'M','T'},-1}, {{'M','V'},1}, {{'M','W'},-1}, {{'M','Y'},-1},
	{{'N','A'},-2}, {{'N','C'},-3}, {{'N','D'},1}, {{'N','E'},0}, {{'N','F'},-3},
	{{'N','G'},0}, {{'N','H'},1}, {{'N','I'},-3}, {{'N','K'},0}, {{'N','L'},-3},
	{{'N','M'},-2}, {{'N','N'},6}, {{'N','P'},-2}, {{'N','Q'},0}, {{'N','R'},0},
	{{'N','S'},1}, {{'N','T'},0}, {{'N','V'},-3}, {{'N','W'},-4}, {{'N','Y'},-2},
	{{'P','A'},-1}, {{'P','C'},-3}, {{'P','D'},-1}, {{'P','E'},-1}, {{'P','F'},-4},
	{{'P','G'},-2}, {{'P','H'},-2}, {{'P','I'},-3}, {{'P','K'},-1}, {{'P','L'},-3},
	{{'P','M'},-2}, {{'P','N'},-2}, {{'P','P'},7}, {{'P','Q'},-1}, {{'P','R'},-2},
	{{'P','S'},-1}, {{'P','T'},-1}, {{'P','V'},-2}, {{'P','W'},-4}, {{'P','Y'},-3},
	{{'Q','A'},-1}, {{'Q','C'},-3}, {{'Q','D'},0}, {{'Q','E'},2}, {{'Q','F'},-3},
	{{'Q','G'},-2}, {{'Q','H'},0}, {{'Q','I'},-3}, {{'Q','K'},1}, {{'Q','L'},-2},
	{{'Q','M'},0}, {{'Q','N'},0}, {{'Q','P'},-1}, {{'Q','Q'},5}, {{'Q','R'},1},
	{{'Q','S'},0}, {{'Q','T'},-1}, {{'Q','V'},-2}, {{'Q','W'},-2}, {{'Q','Y'},-1},
	{{'R','A'},-1}, {{'R','C'},-3}, {{'R','D'},-2}, {{'R','E'},0}, {{'R','F'},-3},
	{{'R','G'},-2}, {{'R','H'},0}, {{'R','I'},-3}, {{'R','K'},2}, {{'R','L'},-2},
	{{'R','M'},-1}, {{'R','N'},0}, {{'R','P'},-2}, {{'R','Q'},1}, {{'R','R'},5},
	{{'R','S'},-1}, {{'R','T'},-1}, {{'R','V'},-3}, {{'R','W'},-3}, {{'R','Y'},-2},
	{{'S','A'},1}, {{'S','C'},-1}, {{'S','D'},0}, {{'S','E'},0}, {{'S','F'},-2},
	{{'S','G'},0}, {{'S','H'},-1}, {{'S','I'},-2}, {{'S','K'},0}, {{'S','L'},-2},
	{{'S','M'},-1}, {{'S','N'},1}, {{'S','P'},-1}, {{'S','Q'},0}, {{'S','R'},-1},
	{{'S','S'},4}, {{'S','T'},1}, {{'S','V'},-2}, {{'S','W'},-3}, {{'S','Y'},-2},
	{{'T','A'},0}, {{'T','C'},-1}, {{'T','D'},-1}, {{'T','E'},-1}, {{'T','F'},-2},
	{{'T','G'},-2}, {{'T','H'},-2}, {{'T','I'},-1}, {{'T','K'},-1}, {{'T','L'},-1},
	{{'T','M'},-1}, {{'T','N'},0}, {{'T','P'},-1}, {{'T','Q'},-1}, {{'T','R'},-1},
	{{'T','S'},1}, {{'T','T'},5}, {{'T','V'},0}, {{'T','W'},-2}, {{'T','Y'},-2},
	{{'V','A'},0}, {{'V','C'},-1}, {{'V','D'},-3}, {{'V','E'},-2}, {{'V','F'},-1},
	{{'V','G'},-3}, {{'V','H'},-3}, {{'V','I'},3}, {{'V','K'},-2}, {{'V','L'},1},
	{{'V','M'},1}, {{'V','N'},-3}, {{'V','P'},-2}, {{'V','Q'},-2}, {{'V','R'},-3},
	{{'V','S'},-2}, {{'V','T'},0}, {{'V','V'},4}, {{'V','W'},-3}, {{'V','Y'},-1},
	{{'W','A'},-3}, {{'W','C'},-2}, {{'W','D'},-4}, {{'W','E'},-3}, {{'W','F'},1},
	{{'W','G'},-2}, {{'W','H'},-2}, {{'W','I'},-3}, {{'W','K'},-3}, {{'W','L'},-2},
	{{'W','M'},-1}, {{'W','N'},-4}, {{'W','P'},-4}, {{'W','Q'},-2}, {{'W','R'},-3},
	{{'W','S'},-3}, {{'W','T'},-2}, {{'W','V'},-3}, {{'W','W'},11}, {{'W','Y'},2},
	{{'Y','A'},-2}, {{'Y','C'},-2}, {{'Y','D'},-3}, {{'Y','E'},-2}, {{'Y','F'},3},
	{{'Y','G'},-3}, {{'Y','H'},2}, {{'Y','I'},-1}, {{'Y','K'},-2}, {{'Y','L'},-1},
	{{'Y','M'},-1}, {{'Y','N'},-2}, {{'Y','P'},-3}, {{'Y','Q'},-1}, {{'Y','R'},-2},
	{{'Y','S'},-2}, {{'Y','T'},-2}, {{'Y','V'},-1}, {{'Y','W'},2}, {{'Y','Y'},7},
};

SeqAlignment SeqAlignment::do_Needleman(const std::string & seqa,const std::string & seqb,
		const SMatrix & smatrix,double open, double ext){
	SeqAlignment Al;
	Al.seqa_=seqa;
	Al.seqb_=seqb;
	std::vector<std::vector<double>> sij;
	std::vector<std::vector<std::pair<int,int>>> trace;
	sij.resize(seqa.size()+1);
	trace.resize(seqa.size()+1);
	for(auto &row:sij)row.resize(seqb.size()+1);
	for(auto &row:trace) row.resize(seqb.size()+1);
	std::vector<double> & s0=sij.at(0);
	std::vector<std::pair<int,int>> & t0=trace.at(0);
	s0[0]=0;
	t0[0]=std::make_pair(-1,-1);
	for(int j=1;j<seqb.size()+1; ++j){
		s0[j]=open + (double) (j-1) * ext;
		t0[j]=std::make_pair(0,0);
	}
	for(int i=1; i<seqa.size()+1; ++i) {
		sij[i][0]=open+(double) (i-1)*ext;
		trace[i][0]=std::make_pair(0,0);
	}
	for(int i=1; i<seqa.size()+1;++i) {
		std::vector<double> &si=sij.at(i);
		std::vector<std::pair<int,int>> &ti=trace.at(i);
		for (int j=1;j<seqb.size()+1;++j) {
			auto eleiter=smatrix.find(std::make_pair(seqa[i-1],seqb[j-1]));
			double sab;
			if(eleiter == smatrix.end()) sab=ext;
			else sab=eleiter->second;
			double smax=sij[i-1][j-1]+sab;
			int tu=i-1;
			int tl=j-1;
			for(int k=0;k<i; ++k) {
				double s_g=sij[k][j] + open + (double) (i-k-1) *ext;
				if(s_g>smax) {
					smax=s_g;
					tu=k;
					tl=j;
				}
			}
			for(int k=0; k<j; ++k) {
				double s_g=si[k] + open +double(j-k-1)*ext;
				if(s_g>smax){
					smax=s_g;
					tu=i;
					tl=k;
				}
			}
			si[j]=smax;
			ti[j]=std::make_pair(tu,tl);
		}
	}
	double totalscore=sij[seqa.size()][seqb.size()];
	std::vector<std::pair<int,int>> path;
	int ic=seqa.size();
	int jc=seqb.size();
	path.push_back(std::make_pair(ic,jc));
	bool cont_trace=true;
	while(cont_trace) {
		std::pair<int,int> & ij_new=trace[ic][jc];
		ic=ij_new.first;
		jc=ij_new.second;
		path.push_back(std::make_pair(ic,jc));
		cont_trace=(ic >0 || jc >0);
	}
	Al.b2a_.resize(seqa.size(),-1);
	Al.a2b_.resize(seqb.size(),-1);
	Al.aligned_.clear();
	for(int p=path.size()-1; p>0; --p) {
		int ip=path[p].first;
		int jp=path[p].second;
		int i=path[p-1].first;
		int j=path[p-1].second;
		if(i-ip==1 && j-jp == 1){
			Al.b2a_[i-1]=j-1;
			Al.a2b_[j-1]=i-1;
			Al.aligned_.push_back(std::make_pair(i-1,j-1));
		} else if(i-ip > 0) {
			assert(j==jp);
			for(int k=ip+1;k<=i;++k) {
				Al.aligned_.push_back(std::make_pair(k-1,-1));
			}
		} else {
			assert(i==ip);
			for(int k=jp+1; k<=j; ++k){
				Al.aligned_.push_back(std::make_pair(-1,k-1));
			}
		}
	}
	return Al;
}

