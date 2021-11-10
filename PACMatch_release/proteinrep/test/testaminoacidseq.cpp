/*
 * testaminoacidseq.cpp
 *
 *  Created on: 2017年6月17日
 *      Author: hyliu
 */

#include "proteinrep/aminoacidseq.h"
#include <iostream>
using namespace NSPproteinrep;
int main(int argc, char **argv){
	std::string seq1{"YTAEQLMCTVYWGLTVEQGAL"};
	std::string seq2{"ADQXLTVYWEGLTVEQMELQQ"};
	SeqAlignment al=SeqAlignment::do_Needleman(seq1,seq2,SeqAlignment::blosum62matrix,0,0);
	std::cout << al.alignedseqa()<<std::endl;
	std::cout << al.alignedseqb()<<std::endl;
	std::cout <<"identity: " <<al.identity() <<std::endl;
	for(int i=0;i<seq1.size();++i)
		std::cout <<' '<<al.alignedposi_b(i);
	std::cout <<std::endl;
	for(int i=0;i<seq2.size();++i)
		std::cout <<' '<<al.alignedposi_a(i);
	std::cout <<std::endl;

}

