/*
 * denovopocketdesign.cpp
 *
 *  Created on: 2020年7月14日
 *      Author: yxchen
 */
/* new version
 * 1. read all bf_contained template.sdf (not protein.pdb)
 * and decompose them into common_motif(target_atms), vector<OBMol>.
 * 2. from largest common_motif to smallest, make ssite. For the same common_motif,
 * pick top perc. ssite. This process until all atoms_needs meet or all templates screened.
 * 3. build pocket from ssites.
*/

#include "buildpocket.h"
#include "bfdecoupling.h"
#include <ctime>
#include <time.h>
#include <dataio/datapaths.h>
using namespace myobcode;
using namespace subsitedesign;
using namespace OpenBabel;
using namespace NSPproteinrep;
using namespace buildpocket;
using namespace std;
using namespace NSPdataio;
#define DIST2CLASH 25.0

int main(int argc, char **argv)
{
	ofstream ofstr("design_record.txt");
	clock_t startT, Tnow;
	startT = clock();
// read par
    ofstr << "Reading parameters ..." << endl;
	std::string par(argv[1]);
	std::ifstream ifpar;
	ifpar.open(par.c_str());
	if (!ifpar.good())
	{
		ofstr << "no par" << std::endl;
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
// some parameters
	int seed = std::stoi(getp(parmap, "seed"));
	srand((unsigned)seed);
	int needs = std::stoi(getp(parmap, "bf_needs"));
	double perc = std::stod(getp(parmap, "percentage")); // all tmplts should have simi < perc.
	double th = std::stod(getp(parmap, "threshold"));
	double wdmldpl = std::stod(getp(parmap, "weight_dmldpl"));

// read target_bf.txt to make bf_names
	std::string target_bf = getp(parmap, "target_bf");
	std::ifstream ifs;
	ifs.open(target_bf.c_str());
	if (!ifs.good())
	{
		ofstr << "sdf failure" << std::endl;
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

// read & prepare target
    ofstr << "Decomposing target ligand ..." << endl;
	const std::map<std::string, BasicFragment> &map =
				BasicFragment::getmap();
	OBConversion obconversion;
	OBMol targetmol;
	OBMol templatemol;
	obconversion.SetInFormat("sdf");
	std::string targetsdf = getp(parmap, "target.sdf");
	obconversion.ReadFile(&targetmol, targetsdf);
	//determine atom types of target ligand
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

	std::map<std::string, std::vector<std::shared_ptr<FMMatches>>> fragtargetmatches;
	for (auto &m : map)
		fragtargetmatches[m.first] = findfmmatches(&(m.second), &targetmol, false);

// build, rank & fetch top_perc. template OBMols by similar
	// bf_names -> topperc.simi_names
	std::vector<std::vector<std::string>> bf_names; // [bflists]][tmplt_names]
	for (int b = 0; b < bflists.size(); b++)
	{
		std::string bf_fn = getenvpath("DEPACT_DATAPATH")+"bf_"+bflists[b]+".txt";
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
	std::vector<OBMol> tmplt_mols;
	topsimitmplts(tmplt_mols, targetmol, bf_names, needs, perc);
	Tnow = clock();
	ofstr << "Time: " << (double)(Tnow-startT)/(3600*CLOCKS_PER_SEC) << " hours"
			<< " Fetch top_similar tmplts " << tmplt_mols.size() << std::endl;

	// decompose into ssas.
	std::vector<TmpltSSAs> tmpltssas;
	for (auto tm : tmplt_mols)
	{
		std::map<std::string, std::vector<std::shared_ptr<FMMatches>>> fragtmpltmatches;
		for (auto &m : map)
		fragtmpltmatches[m.first] = findfmmatches(&(m.second), &tm, true);
		//determine atomtypes of template
		std::vector<std::vector<std::string>> tmpltatomtypes_all=findatomtypes(&tm);
		std::vector<std::string> tmpltatomtypes;
		bool knownatoms = true; // require all atoms in this tm should be known.
		for(auto &tat:tmpltatomtypes_all)
		{
			if(tat.empty())
			{
				knownatoms = false;//tmpltatomtypes.push_back("unspecified");
				break;
			}
			else tmpltatomtypes.push_back(tat[0]);
		}
		if (!knownatoms)
		{
			std::cout << tm.GetTitle() << " is passed for having unspecified atom(s)." << std::endl;
			continue;
		}
		//find substructure alignments
		std::vector<std::shared_ptr<SubstrAlignment>> ssas;
		ssas = findssas(fragtargetmatches, fragtmpltmatches);
		if (ssas.size() == 0) continue;
		//wrap up the substruature alignments and other data
		auto tmpltssa = maketmpltssas(ssas,tmpltatomtypes);
		tmpltssas.push_back(*tmpltssa);
	}
	Tnow = clock();
	ofstr << "Time: " << (double)(Tnow-startT)/(3600*CLOCKS_PER_SEC) << " hours"
			<< " Finsih tmplt_alignments." << std::endl;

	//make ssites
	int idx = 0;
	std::vector<std::shared_ptr<Subsite>> pre_ssites; // every ssa-> one ssite
	std::map<int, std::vector<std::string>> covernum_codes;
	std::map<std::string, std::map<double, std::vector<int>>> pcode_ids;
	int tgt_anum = tgt.conformer.atomlist.size();
	for (int i = 0; i < tmpltssas.size(); i++)
	{
		auto &ssa = tmpltssas[i];
		if (i > 0 && i % 1000 == 0)
			ofstr << "makessites: " << i <<std::endl;
//		makessites(ssa, tgt, wdmldpl, pre_ssites, idx, covernum_codes, code_ids);
		// copy makessites because pre_ssites' TargetStruct* is wrong when using subfunction.
		std::vector<std::string> parsedname = TmpltSSAs::parsetmpltname(
				ssa.tmpltname);
		int modelno=std::atoi(parsedname[TmpltSSAs::MODELNO].c_str());
		char chainid=parsedname[TmpltSSAs::CHAINID][0];
		int residuenumber=std::atoi(parsedname[TmpltSSAs::RESIDUENO].c_str());
		auto reader=readpdbmodel(parsedname[TmpltSSAs::PDBID],1,modelno);
		if (reader == nullptr)
			std::cout << "find an empty: " << parsedname[TmpltSSAs::PDBID] << std::endl;
		else
		{
			MapPdbKeyInt & mki = *(reader->mappdbkeyint());
			int ichain = mki.chainNumber(chainid);
			//reskey is pair of resiude sequence and insertion code from PDB record
			int iresidue = mki.posiNumber(std::pair<int,char>(residuenumber,' '), chainid);
			if (iresidue != -10)
			{
				AAConformersInModel cim;
				cim.getconformers(*reader);
				auto contacts = findcontacts(&cim,ichain,iresidue);
				for (auto &aln : ssa.alignments)
				{
					std::shared_ptr<Subsite> ss(new Subsite(tgt,*contacts,aln));
					if (ss->keyresidues().empty()) continue;
					int codenum = 0;
					string code(tgt_anum, '0'); // empty_atm_code
					for (auto ap : aln.alignedpairs)
					{
						code[ap.first] = '1'; // ap.first has already -1 during ssas->tmpltssas
						codenum--;
					}
					double tscore = ss->totalscore(wdmldpl);
					if (pcode_ids.count(code) > 0)
					{
						if (pcode_ids[code].count(tscore) > 0)
						{
							// because same code & same score ~= totally the same
							continue;
//							pcode_ids[code][tscore].push_back(idx);
//							std::cout << "find same code & score ssites in " << parsedname[TmpltSSAs::PDBID] << std::endl;
						}
						else //if (pcode_ids[code].count(tscore) == 0)
						{
							vector<int> idxs {idx};
							pcode_ids[code].insert(make_pair(tscore, idxs));
						}
					}
					else
					{
						std::map<double, std::vector<int>> score_idxs;
						std::vector<int> idxs {idx};
						score_idxs.insert(make_pair(tscore, idxs));
						pcode_ids.insert(make_pair(code, score_idxs));
					}
					if (covernum_codes.count(codenum) > 0)
						covernum_codes[codenum].push_back(code);
					else
					{
						vector<string> codes {code};
						covernum_codes.insert(make_pair(codenum, codes));
					}
					pre_ssites.push_back(ss);
					idx++;
				}
			}
			else
				std::cout << "wrong resnum " << residuenumber << " in " << chainid << std::endl;
		} // makessites.copy
	} // for every tmpltssas -> pre_ssites

	Tnow = clock();
	ofstr << "Time: " << (double)(Tnow-startT)/(3600*CLOCKS_PER_SEC) << " hours"
			<< " pre_ssites.size() " << pre_ssites.size() << std::endl;

	// pre_ssites -> unredundent_ssites: Atom_different, otherwise RMSD > 0.3
//	std::vector<std::shared_ptr<Subsite>> ussites; // top_score_ssites.
//	std::map<std::string, std::map<double, std::vector<int>>> ucode_ids;
//	unredundent(pre_ssites, ussites, pcode_ids, ucode_ids, 0.3);
	// unredundent_ssites -> ssites.
	std::vector<std::shared_ptr<Subsite>> ssites; // top_score_ssites.
	std::map<std::string, std::map<double, std::vector<int>>> code_ids;
	topssites(pre_ssites, ssites, pcode_ids, code_ids, th);

	Tnow = clock();
	ofstr << "Time: " << (double)(Tnow-startT)/(3600*CLOCKS_PER_SEC) << " hours"
			<< " Finish building ssites: " << ssites.size() << std::endl;

// fetch special contacts, which is forced to add into the pocket, such as coenzyme
	std::vector<OBMol> scmols;
	auto iter_par = parmap.find("special_contacts");
	if (iter_par == parmap.end())
	{
		ofstr << "need special_contacts in par file" << std::endl;
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
	std::vector<int> s2f(ssites.size(), -1); // ssites_id TO fssites_id; should have a test.
	int sf = 0; // sid
	int f = 0; // fid
	std::vector<std::shared_ptr<Subsite>> fssites;
	if (sccrds.size() > 0)
	{
		for (auto site : ssites)
		{
			if (compactsc(sccrds, *site))
			{
				fssites.push_back(site);
				s2f[sf++] = f++;
			}
			else
				sf++;
		}
	}
	else
		for (auto site : ssites)
		{
			fssites.push_back(site);
			s2f[sf++] = f++;
		}

	std::map<std::string, std::map<double, std::vector<int>>> fcode_ids; // code_ids for fssites
	for (auto ita = code_ids.begin(); ita != code_ids.end(); ita++)
	{
		std::map<double, std::vector<int>> fci;
		for (auto itb = ita->second.begin(); itb != ita->second.end(); itb++)
		{
			std::vector<int> fi;
			for (auto i : itb->second)
				if (s2f[i] != -1)
					fi.push_back(s2f[i]);
			if (fi.size() > 0)
				fci.insert(std::make_pair(itb->first, fi));
		} // build fci
		if (!fci.empty())
			fcode_ids.insert(std::make_pair(ita->first, fci));
	} // build fcode_ids
	Tnow = clock();
	ofstr << "Time: " << (double)(Tnow-startT)/(3600*CLOCKS_PER_SEC) << " hours"
			<< " Final ssites compatible: " << fssites.size() << std::endl;
	if (fssites.size() == 0) exit(1);

//check fssites
	iter_par = parmap.find("output_fssites");
	if (iter_par != parmap.end())
		if (iter_par->second == "YES")
		{
			obconversion.SetOutFormat("pdb");
			int ssidx=0;
			for (auto st : fssites)
			{
				std::string outfile="fssites_"+std::to_string(ssidx++)+".pdb";
				st->writepdb(outfile);
			}

			// test: check fcode_ids
			ofstream offc("fcode_ids.txt");
			offc << "code score fssite_id" << std::endl;
			for (auto itf = fcode_ids.begin(); itf != fcode_ids.end(); itf++)
				for (auto itfc = itf->second.begin(); itfc != itf->second.end(); itfc++)
					for (auto k : itfc->second)
						offc << itf->first << " " << itfc->first << " " << to_string(k) << std::endl;

		}
		else if (iter_par->second != "NO")
		{
			ofstr << "need output_ssites in par file to be YES or NO" << std::endl;
			exit(1);
		}


//ver.1 metropolis fetch fssites to make pocket
	int p_num = std::stoi(getp(parmap, "pocket_num"));
	int no_change_num;
	iter_par = parmap.find("no_change_num");
	if (iter_par != parmap.end())
		if (iter_par->second == "auto") no_change_num = floor(fssites.size() * 0.4); // smaller no_change_num -> more variety
		else no_change_num = std::stoi(iter_par->second);
	else
	{
		ofstr << "need no_change_num in par file" << std::endl;
		exit(1);
	}
	int tot_run_num;
	iter_par = parmap.find("tot_run_num");
	if (iter_par != parmap.end())
		if (iter_par->second == "auto") tot_run_num = 5*fssites.size();
		else tot_run_num = std::stoi(iter_par->second);
	else
	{
		ofstr << "need tot_run_num in par file" << std::endl;
		exit(1);
	}
	std::map<std::string, std::vector<int>> pockets_contacts;
	std::vector<std::shared_ptr<Subsite>> final_pockets;
	std::ofstream ofs("pocket_score.txt");
	std::map<int, double> pscores; // pid, score
//	std::map<int, int> pcovers; // pid, atm_cover
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
			else if (1/exp(0.466*(score-sums[minp])) > rand()/(RAND_MAX+0.0)) // 1/e^(0.47*1.5) ~ 0.5
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
//		std::string pd = "pocket_details.txt";
//		std::vector<double> tscores = pocket->detailtotalscore(pd, wdmldpl);
		double ts = pocket->totalscore(wdmldpl);
		ofs << "Pocket_" << p << " total score: " << ts << std::endl;

		ofs << " mmber: ";
		for (auto &m : min_mmber)
			ofs << m << " ";
		ofs << std::endl;
	// ver.1 end

/*
// ver.2 code_dependent_MC
	int p_num = std::stoi(getp(parmap, "pocket_num"));
	int no_change_num; // code_no_change num.
	iter_par = parmap.find("no_change_num");
	if (iter_par != parmap.end())
		if (iter_par->second == "auto") no_change_num = floor(fssites.size() * 0.4); // smaller no_change_num -> more variety
		else no_change_num = std::stoi(iter_par->second);
	else
	{
		ofstr << "need no_change_num in par file" << std::endl;
		exit(1);
	}
	int tot_run_num; // total_runs
	iter_par = parmap.find("tot_run_num");
	if (iter_par != parmap.end())
		if (iter_par->second == "auto") tot_run_num = 5*fssites.size();
		else tot_run_num = std::stoi(iter_par->second);
	else
	{
		ofstr << "need tot_run_num in par file" << std::endl;
		exit(1);
	}
	std::map<std::string, std::vector<int>> pockets_contacts;
	std::vector<std::shared_ptr<Subsite>> final_pockets;
	std::ofstream ofs("pocket_score.txt");
	// test
	std::ofstream oftest("process_test.txt");

	std::map<int, double> pscores; // pid, score
//	std::map<int, int> pcovers; // pid, atm_cover

	int tot_pick = fssites.size();
	for (int p = 0; p < p_num; p++)
	{

		oftest << "Pocket " << p << std::endl;

		double score = 0.0;
		std::vector<int> old_mmber;
		double min_score = 0.0;
		std::vector<int> min_mmber;
		int run = 0;
		int no_change = 0;
		int tp = 0;
		while (run < tot_run_num)
		{
			string t_code(tgt_anum, '1'); // total_code.
			oftest << "run " << to_string(run) << ": " << t_code << std::endl;

			while (no_change < no_change_num)
			{
				// rand_id
				vector<string> max_score_codes = maxscorecodes(fcode_ids, t_code);
				int msc = rand() % max_score_codes.size(); // [0, size()-1]
				auto fci = fcode_ids[max_score_codes[msc]];
				int id = unirandid(fci);

				oftest << "random pick " << id << " " << msc << " " << max_score_codes[msc] << std::endl;

				tp++;
				if (tp > tot_pick) break;
				// renew min_mmber
				std::vector<std::vector<int>> mmbers;
				std::vector<std::shared_ptr<Subsite>> new_subss;

				for (auto m : old_mmber)
				{
					new_subss.push_back(fssites[m]);
					oftest << m << " "; // test
				}
				oftest << std::endl;

				new_subss.push_back(fssites[id]);
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
				else if (exp(0.47*(score-sums[minp])) > rand()/(RAND_MAX+0.0)) // 1/e^(0.47*1.5) ~ 0.5
				{
					std::vector<int> new_mmber;
					for (auto m : mmbers[minp])
					{
						assert(m <= old_mmber.size());
						if (m == old_mmber.size())
							new_mmber.push_back(id);
						else
							new_mmber.push_back(old_mmber[m]);
					}
					old_mmber.clear();
					for (auto m : new_mmber)
						old_mmber.push_back(m);
					score = sums[minp];
					no_change = 0;
					for (int i = 0; i < t_code.size(); i++)
						if (t_code[i] == '1' && max_score_codes[msc][i] == '1')
							t_code[i] = '0';
				}
				else
					no_change++;
				//renew min_mmber
				if (score < min_score)
				{
					min_score = score;
					min_mmber.clear();
					for (auto m : old_mmber)
						min_mmber.push_back(m);
				}
				bool tc_finish = true;
				for (auto tc : t_code)
					if (tc == '1')
					{
						tc_finish = false;
						break;
					}
				if (tc_finish) break;
			} // no_change
			if (tp > tot_pick) break;
			run++;
		} // run

		std::shared_ptr<Subsite> pocket = fssites.at(min_mmber[0]);
		for (int i = 1; i < min_mmber.size(); ++i)
			pocket = std::shared_ptr <Subsite> (new Subsite(*pocket, *(fssites.at(min_mmber[i]))));

		final_pockets.push_back(pocket);
		std::string file = "Pocket_"+std::to_string(p)+".pdb";
		pocket->writepdb(file);
		std::string pd = "pocket_details.txt";
		std::vector<double> tscores = pocket->detailtotalscore(pd, wdmldpl);
		ofs << "Pocket_" << p << " total score: " << tscores[0] << std::endl;
		ofs << " mmber: ";
		for (auto &m : min_mmber)
			ofs << m << " ";
		ofs << std::endl;
	//ver.2 end
*/
// record score for every pocket
		pscores.insert(make_pair(p, ts));
	} // every pocket
// output pockets' evaluation
	ofs << "Pockets Evaluation:" << std::endl;
	ofs << "Pocket_ID  Score" << std::endl;
	for (int p = 0; p < p_num; p++)
		ofs << p << " " << pscores[p] << std::endl;
//		ofs << p << " " << pscores[p] << " " << pcovers[p] << std::endl;

	Tnow = clock();
	ofstr << "Time: " << (double)(Tnow-startT)/(3600*CLOCKS_PER_SEC) << " hours"
			<< " Finish testoverall.cpp" << std::endl;
}



