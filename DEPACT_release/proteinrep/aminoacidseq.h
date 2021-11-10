/*
 * aminoacidseq.h
 *
 *  Created on: 2017年6月16日
 *      Author: hyliu
 */

#ifndef PROTEINREP_AMINOACIDSEQ_H_
#define PROTEINREP_AMINOACIDSEQ_H_
#include <string>
#include <map>
#include <vector>

namespace NSPproteinrep {
class AminoAcidSeq {
public:
	typedef std::vector<std::string> NameSeq;
	typedef std::string CodeSeq;
	static std::map<char,std::string> CodeName;
	static std::map<std::string,char> NameCode;
	static std::string code2name(char c){
		if(CodeName.find(c) != CodeName.end()) return CodeName.at(c);
		return "ANY";
	}
	static char name2code(const std::string & name) {
		if(NameCode.find(name) != NameCode.end())return NameCode[name];
		return 'X';
	}
	static std::string name2code(const std::vector<std::string> &nameseq);
	static std::vector<std::string> code2name(const std::string &codeseq);
};

class SeqAlignment {
public:
	typedef std::map<std::pair<char,char>,double> SMatrix;
	static SMatrix blosum62matrix;
	static SeqAlignment do_Needleman(const std::string &seqa, const std::string &seqb,
		const SMatrix & matrix=blosum62matrix,double open=-10.0,double ext=-0.5);
	std::string alignedseqa() const{
		std::string res;
		for(auto & p:aligned_){
			if(p.first<0) res.push_back('-');
			else res.push_back(seqa_[p.first]);
		}
		return res;
	}
	std::string alignedseqb() const{
		std::string res;
		for(auto&p:aligned_){
			if(p.second<0) res.push_back('-');
			else res.push_back(seqb_[p.second]);
		}
		return res;
	}
	double identity() const {
		int nid=0;
		for(auto & p:aligned_){
			if(p.first<0 || p.second<0) continue;
			if(seqa_[p.first] == seqb_[p.second]) ++nid;
		}
		int maxlen=seqa_.size()>seqb_.size()? seqa_.size():seqb_.size();
		return (double)nid/(double) maxlen;
	}
	double coveratio_a() const {
		int a_cover=0;
		for(auto pb:b2a_) if(pb>=0) ++a_cover;
		return (double) a_cover/ (double) seqa_.size();
	}
	double coveratio_b() const {
		int b_cover=0;
		for(auto pa:a2b_) if(pa>=0) ++b_cover;
		return (double) b_cover/(double) seqb_.size();
	}
	int alignedlength() const {return aligned_.size();}
	int alignedposi_b(int posia) const{
		return b2a_.at(posia);
	}
	int alignedposi_a(int posib) const{
		return a2b_.at(posib);
	}
private:
	std::string seqa_;
	std::string seqb_;
	std::vector<int> b2a_;
	std::vector<int> a2b_;
	std::vector<std::pair<int,int>> aligned_;
};
}


#endif /* PROTEINREP_AMINOACIDSEQ_H_ */
