/*
 * testoverall.cpp
 *
 *  Created on: 2019年2月21日
 *      Author: yxchen
 */

// the correct & new version is denovopocketdesign.cpp

#include "buildpocket.h"
using namespace myobcode;
using namespace subsitedesign;
using namespace OpenBabel;
using namespace NSPproteinrep;
using namespace buildpocket;
#define DIST2CLASH 25.0

#include "bfdecoupling.h"


int main(int argc, char **argv)
{
// read par
	std::string par(argv[1]);
	std::ifstream ifpar;
	ifpar.open(par.c_str());
	if (!ifpar.good())
	{
		std::cout << "no par" << std::endl;
		exit(1);
	}
	std::map<std::string, std::string> parmap;
	while(true)
	{
		std::string line;
		line.clear();
		getline(ifpar, line);
		if(!ifpar.good()) break;
		if (line.size()==0) continue;
		std::vector<std::string> words;
		std::stringstream input(line);
		std::string word;
		while(input>>word) words.push_back(word);
		parmap.insert(std::pair<std::string, std::string>(words[0], words[2]));
	}
	ifpar.close();

// read target_bf.txt to make bf_names
	std::string target_bf = getp(parmap, "target_bf");
	std::ifstream ifs;
	ifs.open(target_bf.c_str());
	if (!ifs.good())
	{
		std::cout << "sdf failure" << std::endl;
		exit(1);
	}
	int start = 1;
	std::string target_name;
	std::vector<std::string> bflists;
	while(true)
	{
		std::string line;
		line.clear();
		getline(ifs, line);
		if(!ifs.good()) break;
		if (line.size()==0) continue;
		if (start == 1)
		{
			target_name = line;
			start = 0;
		}
		else
			bflists.push_back(line);
	}
	ifs.close();

	std::vector<std::vector<std::string>> bf_names; // [bflists]][tmplt_names]
	for (int b = 0; b < bflists.size(); b++)
	{
		std::string bf_fn = "bf_"+bflists[b]+".txt";
		std::ifstream if_bf;
		if_bf.open(bf_fn.c_str());
		std::vector<std::string> names;
		while(true)
		{
			std::string line;
			line.clear();
			getline(if_bf, line);
			if(!if_bf.good()) break;
			if (line.size()==0) continue;
			names.push_back(line);
		}
		if_bf.close();
		bf_names.push_back(names);
	}

// read target
	const std::map<std::string, BasicFragment> &map =
				BasicFragment::getmap();
	OBConversion obconversion;
	OBMol targetmol;
	OBMol templatemol;
	obconversion.SetInFormat("sdf");
	std::string targetsdf = getp(parmap, "target.sdf");
	obconversion.ReadFile(&targetmol, targetsdf);
// this part will be delete_begin
/*
	std::map<std::string, std::vector<std::shared_ptr<FMMatches>>> fragtargetmatches;
	for (auto &m : map)
		fragtargetmatches[m.first] = findfmmatches(&(m.second), &targetmol,
				false);
	///determine atom types of target ligand
	std::vector<std::vector<std::string>> targetatomtypes_all=findatomtypes(&targetmol);
	std::vector<std::string> targetatomtypes;
	for(auto &tat:targetatomtypes_all)
	{
		if(tat.empty()) targetatomtypes.push_back("unspecified");
		else targetatomtypes.push_back(tat[0]);
	}
	std::ostringstream pdbsstr;
	obconversion.SetOutFormat("pdb");
	obconversion.Write(&targetmol,&pdbsstr);
	std::istringstream isstrpdb(pdbsstr.str());
	TargetStruct tgt(targetmol.GetTitle(),targetatomtypes,isstrpdb);
*/
// this part will be delete_end
// build template OBMols
	int seed = std::stoi(getp(parmap, "seed"));
	srand((unsigned)seed);
	int needs = std::stoi(getp(parmap, "bf_needs"));
	double perc = std::stod(getp(parmap, "percentage"));
	double thread = std::stod(getp(parmap, "threshold"));
	double wdmldpl = std::stod(getp(parmap, "weight_dmldpl"));

/* new version
 * 1. read all bf_contained template.sdf (not protein.pdb)
 * and decompose them by interested_bf.
 * 2. ranked them by intested_bf_num, and build ssites
 * from max_num (some bfs) to min_num (1 bf), until all bf have enough templates.
 * 3. build pocket from ssites.
*/
	// read all tmplts < simi_perc & rank by bf_num
	std::set<std::string> tmplt_names;
	for (auto bfn : bf_names)
		tmplt_names.insert(bfn.begin(), bfn.end());
	std::cout << "Finish reading all tmplt_names: " << tmplt_names.size() << std::endl;
	std::map<int, std::vector<OBMol>> tmpltmols = readalltmplts(tmplt_names, obconversion,
			targetmol, perc);
	std::cout << "Finish reading tmpltmols." << std::endl;
	// fetch ssites from large key to small until meets bf_needs.
//	std::vector<std::shared_ptr<Subsite>> ssites = tmpltssites(tmpltmols, targetmol, needs, thread, wdmldpl);
//	std::cout << "Finish establishing ssites: " << ssites.size() << std::endl;

	// copy tmpltssites
	std::map<std::string, std::vector<std::shared_ptr<FMMatches>>> fragtargetmatches;
	std::vector<int> bfneeds; // total tmplt_needs for target.
	int unfinished = 0;
	for (auto &m : map)
	{
		fragtargetmatches[m.first] = findfmmatches(&(m.second), &targetmol,
				false);
		std::string out = bfdecoupling(&(m.second), &targetmol);
		if (out == "Yes")
		{
			bfneeds.push_back(needs);
			unfinished++;
		}
		else bfneeds.push_back(0);
	}
	///determine atom types of target ligand
	std::vector<std::vector<std::string>> targetatomtypes_all=findatomtypes(&targetmol);
	std::vector<std::string> targetatomtypes;
	for(auto &tat:targetatomtypes_all)
	{
		if(tat.empty()) targetatomtypes.push_back("unspecified");
		else targetatomtypes.push_back(tat[0]);
	}
	std::ostringstream pdbsstr;
	obconversion.SetOutFormat("pdb");
	obconversion.Write(&targetmol,&pdbsstr);
	std::istringstream isstrpdb(pdbsstr.str());
	TargetStruct tgt(targetmol.GetTitle(),targetatomtypes,isstrpdb);

// for tmpltssas establishment
	std::vector<std::shared_ptr<Subsite>> ssites;
	std::vector<TmpltSSAs> tmpltssas;
	for (auto it = tmpltmols.begin(); it != tmpltmols.end(); it++)
	{
		std::cout << "bf_contains: " << std::to_string(it->first) << std::endl;
		std::cout << "OB_contains: " << std::to_string(it->second.size()) << std::endl;
		if (unfinished == 0) break;
		for (auto &tm : it->second)
		{
			if (unfinished == 0) break;
			std::map<std::string, std::vector<std::shared_ptr<FMMatches>>> fragtmpltmatches;
			int n = 0;
			bool find = false;
			for (auto &m : map)
			{
				fragtmpltmatches[m.first] = findfmmatches(&(m.second), &tm, true);
				std::string out = bfdecoupling(&(m.second), &tm);
				if (out == "Yes" && bfneeds[n] != 0)
				{
					bfneeds[n]--;
					find = true;
					if (bfneeds[n] == 0) unfinished--;
				}
				n++;
			}
			if (!find) continue;
			//find substructure alignments
			std::vector<std::shared_ptr<SubstrAlignment>> ssas;
			ssas = findssas(fragtargetmatches, fragtmpltmatches);
			//determine atomtypes of template
			std::vector<std::vector<std::string>> tmpltatomtypes_all=findatomtypes(&tm);
			std::vector<std::string> tmpltatomtypes;
			for(auto &tat:tmpltatomtypes_all)
			{
				if(tat.empty()) tmpltatomtypes.push_back("unspecified");
				else tmpltatomtypes.push_back(tat[0]);
			}
			//wrap up the substruature alignments and other data
			if (ssas.size() == 0) continue;
			auto tmpltssa = maketmpltssas(ssas,tmpltatomtypes);
			tmpltssas.push_back(*tmpltssa);
		}
	}
	std::cout << "Expecting template_pdbs: " << tmpltssas.size() << std::endl;
	//make ssites
	std::map<double, std::vector<int>> test_scores; // score, <idx for test_ssites>
	int idx = 0;
	for (int i = 0; i < tmpltssas.size(); i++)
	{
		auto &ssa = tmpltssas[i];
		if (i % 100 == 0)
			std::cout << "makessites: " << i <<std::endl;
		makessites(ssa, tgt, wdmldpl, idx, test_scores, ssites);
	}
	// copy tmpltssites end

// this part will be delete_begin
/*
// two choice: separate or combine
// separate:
	std::vector<std::vector<OBMol>> tmpgroups = buildtmpmols_sep(bfneeds, bf_names,
			obconversion, targetmol, perc);
	std::vector<std::shared_ptr<Subsite>> ssites;
	for (auto tg : tmpgroups)
	{
		std::vector<std::shared_ptr<Subsite>> pre_ssites;
		std::vector<TmpltSSAs> tmpltssas;
	//establish tmpltssas
		for (auto &templatemol : tg)
		{
			std::map<std::string, std::vector<std::shared_ptr<FMMatches>>> fragtmpltmatches;
			for (auto &m : map)
				fragtmpltmatches[m.first] = findfmmatches(&(m.second), &templatemol, true);

			///find substructure alignments
			std::vector<std::shared_ptr<SubstrAlignment>> ssas;
			ssas = findssas(fragtargetmatches, fragtmpltmatches);
			///determine atomtypes of template
			std::vector<std::vector<std::string>> tmpltatomtypes_all=findatomtypes(&templatemol);
			std::vector<std::string> tmpltatomtypes;
			for(auto &tat:tmpltatomtypes_all)
			{
				if(tat.empty()) tmpltatomtypes.push_back("unspecified");
				else tmpltatomtypes.push_back(tat[0]);
			}
			///wrap up the substruature alignments and other data
			if (ssas.size() == 0) continue;
			auto tmpltssa = maketmpltssas(ssas,tmpltatomtypes);
			tmpltssas.push_back(*tmpltssa);
		}

	//make ssites
		std::map<double, std::vector<int>> test_scores; // score, <idx for test_ssites>
		std::vector<std::shared_ptr<Subsite>> test_ssites;
		int idx = 0;
		for (int i = 0; i < tmpltssas.size(); i++)
		{
			auto &ssa = tmpltssas[i];
			if (i % 100 == 0)
				std::cout << "makessites: " << i <<std::endl;
			makessites(ssa, tgt, wdmldpl, idx, test_scores, test_ssites);
		}
		// rank and get ssites from test_ssites.
		rankformssites(pre_ssites, thread, test_scores, test_ssites);
		for (auto &ps : pre_ssites)
			ssites.push_back(ps);
	}
*/
// this part will be delete_end
	// test
//	obconversion.SetOutFormat("pdb");
//	int sidx=0;
//	for (auto &s : ssites)
//	{
//		std::cout << "writing ssites " << sidx << std::endl;;
//		std::string outfile="ssites_"+std::to_string(sidx++)+".pdb";
//		s->writepdb(outfile);
//	}

// combine:
/*
	std::vector<OBMol> tmpmols = buildtmpmols(bfneeds, bf_names, obconversion, targetmol, perc);
	std::cout << tmpmols.size() << " templates will be analyzed" << std::endl;
	std::vector<TmpltSSAs> tmpltssas;
//establish tmpltssas
	for (auto &templatemol : tmpmols)
	{
		std::map<std::string, std::vector<std::shared_ptr<FMMatches>>> fragtmpltmatches;
		for (auto &m : map)
			fragtmpltmatches[m.first] = findfmmatches(&(m.second), &templatemol, true);

		///find substructure alignments
		std::vector<std::shared_ptr<SubstrAlignment>> ssas;
		ssas = findssas(fragtargetmatches, fragtmpltmatches);
		///determine atomtypes of template
		std::vector<std::vector<std::string>> tmpltatomtypes_all=findatomtypes(&templatemol);
		std::vector<std::string> tmpltatomtypes;
		for(auto &tat:tmpltatomtypes_all)
		{
			if(tat.empty()) tmpltatomtypes.push_back("unspecified");
			else tmpltatomtypes.push_back(tat[0]);
		}
		///wrap up the substruature alignments and other data
		if (ssas.size() == 0) continue;
		auto tmpltssa = maketmpltssas(ssas,tmpltatomtypes);
		tmpltssas.push_back(*tmpltssa);
	}

//make ssites
	int ss_idx = 0;
	std::vector<std::shared_ptr<Subsite>> ssites;

	std::map<double, std::vector<int>> test_scores; // score, <idx for test_ssites>
	std::vector<std::shared_ptr<Subsite>> test_ssites;
	int idx = 0;
	for (int i = 0; i < tmpltssas.size(); i++)
	{
		auto &ssa = tmpltssas[i];
		if (i % 100 == 0)
			std::cout << "makessites: " << i <<std::endl;
		makessites(ssa, tgt, wdmldpl, idx, test_scores, test_ssites);
	}
	// rank and get ssites from test_ssites.
	rankformssites(ssites, thread, test_scores, test_ssites);
*/

// fetch special contacts, which is forced to add into the pocket, such as coenzyme
	std::vector<OBMol> scmols;
	auto iter_par = parmap.find("special_contacts");
	if (iter_par == parmap.end())
	{
		std::cout << "need special_contacts in par file" << std::endl;
		exit(1);
	}
	else
		if (iter_par->second != "NO") // then there should have special contacts
		{
			ifpar.open(par.c_str());
			while(true)
			{
				std::string line;
				line.clear();
				getline(ifpar, line);
				if(!ifpar.good()) break;
				if (line.size()==0) continue;
				std::vector<std::string> words;
				std::stringstream input(line);
				std::string word;
				while(input>>word) words.push_back(word);
				if (words[0] == "special_contacts")
					for (int i = 1; i < words.size()/2; i++)
					{
						obconversion.SetInFormat(words[i*2].c_str());
						OBMol sc;
						obconversion.ReadFile(&sc, words[i*2+1]);
						scmols.push_back(sc);
					}
			}
			ifpar.close();
		}
	std::vector<NSPgeometry::XYZ> sccrds;
	if (scmols.size() > 0)
		for (auto sc : scmols)
		{
			double *c = sc.GetConformer(0);
			for (int a = 0; a < sc.NumAtoms(); a++)
			{
				NSPgeometry::XYZ cts(c[a*3], c[a*3+1], c[a*3+2]);
				sccrds.push_back(cts);
			}
		}
	//filter for ssites: compatible for special_contacts
	std::vector<std::shared_ptr<Subsite>> fssites;
	if (sccrds.size() > 0)
	{
		for (auto site : ssites)
			if (compactsc(sccrds, *site))
				fssites.push_back(site);
	}
	else
		for (auto site : ssites)
			fssites.push_back(site);
	std::cout << "ssites compatible: " << fssites.size() << std::endl;
	if (fssites.size() == 0) exit(1);

//check ssites
	iter_par = parmap.find("output_fssites");
	if (iter_par != parmap.end())
		if (iter_par->second == "yes")
		{
			obconversion.SetOutFormat("pdb");
			int ssidx=0;
			for (auto st : fssites)
			{
				std::string outfile="fssites_"+std::to_string(ssidx++)+".pdb";
				st->writepdb(outfile);
			}
		}
		else if (iter_par->second != "no")
		{
			std::cout << "need output_ssites in par file to be 'yes' or 'no'" << std::endl;
			exit(1);
		}

//metropolis fetch ssites to make pocket
	int p_num = std::stoi(getp(parmap, "pocket_num"));
	int no_change_num;
	iter_par = parmap.find("no_change_num");
	if (iter_par != parmap.end())
		if (iter_par->second == "auto") no_change_num = floor(fssites.size() * 0.4);
		else no_change_num = std::stoi(iter_par->second);
	else
	{
		std::cout << "need no_change_num in par file" << std::endl;
		exit(1);
	}
	int tot_run_num;
	iter_par = parmap.find("tot_run_num");
	if (iter_par != parmap.end())
		if (iter_par->second == "auto") tot_run_num = fssites.size() * 20;
		else tot_run_num = std::stoi(iter_par->second);
	else
	{
		std::cout << "need tot_run_num in par file" << std::endl;
		exit(1);
	}
	std::map<std::string, std::vector<int>> pockets_contacts;
	std::vector<std::shared_ptr<Subsite>> final_pockets;
	std::ofstream ofs("pocket_score.txt");
	std::map<double, std::vector<int>> pscores;
	for (int p = 0; p < p_num; p++)
	{
		double score = 0.0;
		std::vector<int> old_mmber;
		double min_score = 0.0;
		std::vector<int> min_mmber;
		int no_change = 0;
		int tot_run = 0;
		while (true)
		{
			tot_run++;
			std::vector<std::vector<int>> mmbers;
			std::vector<std::shared_ptr<Subsite>> new_subss;
			for (auto m : old_mmber)
				new_subss.push_back(fssites[m]);
			int new_m = floor(rand()/(RAND_MAX+0.0)*fssites.size());
			new_subss.push_back(fssites[new_m]);
			std::vector<std::shared_ptr<Subsite>> subpockets = combinebycliques_mm(new_subss, mmbers);
			//the best pocket for this subss
			std::vector<double> sums;
			for (auto &ss : subpockets)
				sums.push_back(ss->totalscore(wdmldpl));
			int minp = min_element(sums.begin(), sums.end()) - sums.begin();
            //if new_sub makes clash and raise the score,
			//then we should compare this clash one with the old.
			assert(subpockets.size() <= 2);
			if (subpockets.size() == 2 && score == sums[minp])
				minp = 1-minp;

			if (score == sums[minp])
			{
				no_change++;// as == means the original pocket
			}
			else if (exp(0.7*(score-sums[minp])) > rand()/(RAND_MAX+0.0)) // 1/e^(0.7*1) ~ 0.5
			{
				std::vector<int> new_mmber;
				for (auto m : mmbers[minp])
				{
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
			}
			else
				no_change++;
/*
			//test
			std::cout << score << " ";
			for (auto m : old_mmber)
				std::cout << m << " ";
			std::cout << std::endl;
*/
			//renew min_mmber
			if (score < min_score)
			{
				min_score = score;
				min_mmber.clear();
				for (auto m : old_mmber)
					min_mmber.push_back(m);
			}
			if (no_change == no_change_num || tot_run > tot_run_num)//20; 1000: for 181 fssites
				break;
		}

		std::shared_ptr<Subsite> pocket = fssites.at(min_mmber[0]);
		for (int i = 1; i < min_mmber.size(); ++i)
			pocket = std::shared_ptr <Subsite> (new Subsite(*pocket, *(fssites.at(min_mmber[i]))));

		final_pockets.push_back(pocket);
		std::string file = "Pocket_"+std::to_string(p)+".pdb";
		pocket->writepdb(file);
		std::string pd = "pocket_details.txt";
		std::vector<double> tscores = pocket->detailtotalscore(pd, wdmldpl);
		ofs << "Pocket" << p << " total score: " << tscores[0] << std::endl;

		ofs << "mmber: ";
		for (auto &m : min_mmber)
			ofs << m << " ";
		ofs << std::endl;

// record score for every pocket so as to sort them finally.
		if (pscores.find(min_score) != pscores.end())
			pscores[min_score].push_back(p);
		else
		{
			std::vector<int> ps = {p};
			pscores.insert(std::make_pair(min_score, ps));
		}
	}
// sort pockets by score
	for (auto &ps : pscores)
		for (auto p : ps.second)
			ofs << p << " " << ps.first << std::endl;

	std::cout << "Finish testoverall.cpp" << std::endl;
}



