/*
 * overall_pocket.cpp
 *
 *  Created on: 2019年4月1日
 *      Author: yxchen
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>
#include "subsite.h"
#include "proteinrep/aaconformer.h"

using namespace subsitedesign;
using namespace NSPproteinrep;
// so overall_pocket.cpp can be put into noob

void makessites(TmpltSSAs &ssa, TargetStruct &tgt,
		std::vector<std::shared_ptr<Subsite>> &ssites) {
	std::vector<std::string> parsedname = TmpltSSAs::parsetmpltname(
			ssa.tmpltname);
	int modelno=std::atoi(parsedname[TmpltSSAs::MODELNO].c_str());
	char chainid=parsedname[TmpltSSAs::CHAINID][0];
	int residuenumber=std::atoi(parsedname[TmpltSSAs::RESIDUENO].c_str());
	auto reader=readpdbmodel(parsedname[TmpltSSAs::PDBID],1,modelno);
	if (reader == nullptr)
		std::cout << "find an empty: " << parsedname[TmpltSSAs::PDBID] << std::endl;
	else {
		MapPdbKeyInt & mki = *(reader->mappdbkeyint());
		int ichain = mki.chainNumber(chainid);
		//reskey is pair of resiude sequence and insertion code from PDB record
		int iresidue = mki.posiNumber(std::pair<int,char>(residuenumber,' '), chainid);
		if (iresidue != -10) {
			AAConformersInModel cim;
			cim.getconformers(*reader);
			auto contacts = findcontacts(&cim,ichain,iresidue);
			for (auto &aln : ssa.alignments) {
				std::shared_ptr<Subsite> ss(new Subsite(tgt,*contacts,aln));
				if (ss->keyresidues().empty()) continue;
				ssites.push_back(ss);
			}
		} else
			std::cout << "wrong resnum " << residuenumber << " in " << chainid << std::endl;
	}
}

int main(int argc, char **argv) {
//pocket_num seed
// read ssas.txt
	TargetStruct tgt;
	std::vector<TmpltSSAs> tmpltssas;
	std::ifstream ifs;
	ifs.open("ssas.dat");
	{
		boost::archive::text_iarchive ia(ifs);
		ia >> tgt;
	}
	int ntmplts;
	ifs >> ntmplts;
	for (int i = 0; i < ntmplts; i++){
		TmpltSSAs ssa;
		{
			boost::archive::text_iarchive ia(ifs);
			ia >> ssa;
		}
		tmpltssas.push_back(ssa);
	}
	ifs.close();

//make ssites
	int ss_idx = 0;
	std::vector<std::shared_ptr<Subsite>> ssites;

	for (int i = 0; i < tmpltssas.size(); i++) {
		auto &ssa = tmpltssas[i];
		if (i%100 == 0)
			std::cout << "makessites: " << i <<std::endl;
		makessites(ssa, tgt, ssites);
	}
	std::cout << "ssites.size is: " << ssites.size() << std::endl;

//metropolis fetch ssites to make pocket
	int p_num = std::stoi(argv[1]);
	int seed = std::stoi(argv[2]);
	srand((unsigned)seed);
	std::ofstream ofs("pocket_score.txt", std::ios::app);
	for (int p = 0; p < p_num; p++) {
		double score = 0.0;
		std::vector<int> old_mmber;
		int no_change = 0;
		int tot_run = 0;
		while (true) {
			tot_run++;
			std::cout << tot_run << std::endl;
			std::vector<std::vector<int>> mmbers;
			std::vector<std::shared_ptr<Subsite>> new_subss;
			for (auto m : old_mmber)
				new_subss.push_back(ssites[m]);
			int new_m = floor(rand()/(RAND_MAX+0.0)*ssites.size());
			new_subss.push_back(ssites[new_m]);
			std::vector<std::shared_ptr<Subsite>> pockets = combinebycliques_mm(new_subss, mmbers);
			//the best pocket for this subss
			std::vector<double> sums;
			for (auto &ss : pockets) {
				double sum = 0.0;
				for (auto & pl : ss->keyresiduescores().plscores)
					sum += ss->keyresiduescores().totalplscore(pl.first);
				for (auto & ml : ss->keyresiduescores().mlscores)
					sum += ss->keyresiduescores().totalplscore(ml.first);
				sums.push_back(sum);
			}
			int minp = min_element(sums.begin(), sums.end()) - sums.begin();
            //if new_sub makes clash and raise the score,
			//then we should compare this clash one with the old.
			assert(pockets.size() <= 2);
			if (pockets.size() == 2 && score == sums[minp])
				minp = 1-minp;

			if (score == sums[minp]) {
				no_change++;// as == means the original pocket
			} else if (exp(5*(score-sums[minp])) > rand()/(RAND_MAX+0.0)) {
				std::vector<int> new_mmber;
				for (auto m : mmbers[minp]) {
					assert(m <= old_mmber.size());
					if (m == old_mmber.size())
						new_mmber.push_back(new_m);
					else
						new_mmber.push_back(old_mmber[m]);
				}
				old_mmber.clear();
				for (auto m : new_mmber)
					old_mmber.push_back(m);
				score = sums[minp];
				no_change = 0;
			} else
				no_change++;
			if (no_change == 200 || tot_run > 100000)//for it won't distrub the normal search
				break;
		}
		std::shared_ptr<Subsite> pocket = ssites.at(old_mmber[0]);
		for (int i = 1; i < old_mmber.size(); ++i) {
			pocket = std::shared_ptr < Subsite
					> (new Subsite(*pocket, *(ssites.at(old_mmber[i]))));
		}
		std::string file = "Pocket_"+std::to_string(p)+".pdb";
		pocket->writepdb(file);
		ofs << "Pocket " << std::to_string(p) << " :pl+ml: " << score << std::endl;
		for (auto m : old_mmber)
			ofs << m << " ";
		ofs << std::endl;
	}

}
