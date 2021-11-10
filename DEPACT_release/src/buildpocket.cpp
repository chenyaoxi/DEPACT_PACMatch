/*
 * buildpocket.cpp
 *
 *  Created on: 2019年7月3日
 *      Author: yxchen
 */

#include "buildpocket.h"
#include "bfdecoupling.h"
#include <dataio/datapaths.h>
using namespace myobcode;
using namespace subsitedesign;
using namespace OpenBabel;
using namespace NSPproteinrep;
using namespace buildpocket;
using namespace NSPdataio;
#define DIST2CLASH 25.0

std::map<std::string, std::string> buildpocket::readpar(std::string parname)
{
	std::ifstream ifpar;
	ifpar.open(parname.c_str());
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
	return parmap;
}

std::string buildpocket::getp(std::map<std::string, std::string> parmap, std::string name)
{
	std::string par;
	std::map<std::string, std::string>::iterator iter_par = parmap.find(name);
	if (iter_par != parmap.end())
		par = iter_par->second;
	else
	{
		std::cout << "need " << name << " in par file" << std::endl;
		exit(1);
	}
	return par;
}

bool buildpocket::similarity(OBMol targetmol, OBMol templatemol, double perc)
{
	OBFingerprint* fptype = OBFingerprint::FindType("FP3");
    std::vector<unsigned> fp_tar;
    if(!fptype->GetFingerprint(&targetmol, fp_tar))
    {
    	std::cout << "target fingerprint error" << std::endl;
    	exit(1);
    }
    std::vector<unsigned> fp_tem;
    if(!fptype->GetFingerprint(&templatemol, fp_tem))
    {
    	std::cout << "template fingerprint error" << std::endl;
    	exit(1);
    }
	// calculate simi
    double simi = OBFingerprint::Tanimoto(fp_tar, fp_tem);
	if (simi > perc) // if perc = 1.0, it means we need native sdf as template, so we use >
		return true;
	else
		return false;
}

double buildpocket::similar(OBMol targetmol, OBMol templatemol)
{
	OBFingerprint* fptype = OBFingerprint::FindType("FP3");
    std::vector<unsigned> fp_tar;
    if(!fptype->GetFingerprint(&targetmol, fp_tar))
    {
    	std::cout << "target fingerprint error" << std::endl;
    	exit(1);
    }
    std::vector<unsigned> fp_tem;
    if(!fptype->GetFingerprint(&templatemol, fp_tem))
    {
    	std::cout << "template fingerprint error" << std::endl;
    	exit(1);
    }
	// calculate simi
    double simi = OBFingerprint::Tanimoto(fp_tar, fp_tem);
    return simi;
}

std::vector<OBMol> buildpocket::buildtmpmols(std::vector<int> bfneeds, std::vector<std::vector<std::string>> bf_names,
		OBConversion obconversion, OBMol targetmol, double perc)
{
	std::vector<OBMol> tmpmols;

	std::string t_file = "template.sdf";
	for (int i = 0; i < bfneeds.size(); i++)
	{
		std::cout << "fetch tmplts for bf " << i << std::endl;
		while (bfneeds[i] > 0)
		{
			OBMol templatemol;
			int n = floor(rand()/(RAND_MAX+0.0)*bf_names[i].size());
			std::string tmp_name = bf_names[i][n];
			int find = 0;
			remove("template.sdf");
			std::ifstream ifall;
			ifall.open("all-sdf.sdf");
			while(true)
			{
				std::string line;
				line.clear();
				getline(ifall, line);
				if(!ifall.good()) break;
				if (line.size() == 0) continue;
				if (line == tmp_name) find = 1;
				if (find == 1)
				{
					std::ofstream of_t;
					of_t.open(t_file.c_str(), std::ios::app);
					of_t << line << std::endl;
				}
				if (find == 1 && line == "$$$$") break;
			}
			ifall.close();
			//for tmplt part
			obconversion.ReadFile(&templatemol, t_file);

			if (!similarity(targetmol, templatemol, perc))
			{
				if (tmpmols.size() == 0)
				{
					tmpmols.push_back(templatemol);
					bfneeds[i] -= 1;
				}
				else
				{
					bool find = false;
					for (auto &tmp : tmpmols)
						if (strcmp(tmp.GetTitle(false), templatemol.GetTitle(false)) == 0)
						{ // in fact, strcmp should compaire new.title with bf_names before.
							find = true;
							bfneeds[i] -= 1;
							break;
						}
					if (!find)
					{
						tmpmols.push_back(templatemol);
						bfneeds[i] -= 1;
					}
				}
			} // !similarity
		} // every names
	} // every bf
	return tmpmols;
}

std::vector<std::vector<OBMol>> buildpocket::buildtmpmols_sep(std::vector<int> bfneeds, std::vector<std::vector<std::string>> bf_names,
		OBConversion obconversion, OBMol targetmol, double perc)
{
	std::vector<std::vector<OBMol>> tmgroup;
	std::string t_file = "template.sdf";

	for (int i = 0; i < bfneeds.size(); i++)
	{
		// for every tmpmols
		std::vector<OBMol> tmpmols;
		std::cout << "fetch tmplts for bf " << i << std::endl;
		while (bfneeds[i] > 0)
		{
			OBMol templatemol;
			int n = floor(rand()/(RAND_MAX+0.0)*bf_names[i].size());
			std::string tmp_name = bf_names[i][n];
			int find = 0;
			remove("template.sdf");
			std::ifstream ifall;
			ifall.open("all-sdf.sdf");
			while(true)
			{
				std::string line;
				line.clear();
				getline(ifall, line);
				if(!ifall.good()) break;
				if (line.size() == 0) continue;
				if (line == tmp_name) find = 1;
				if (find == 1)
				{
					std::ofstream of_t;
					of_t.open(t_file.c_str(), std::ios::app);
					of_t << line << std::endl;
				}
				if (find == 1 && line == "$$$$") break;
			}
			ifall.close();
			//for tmplt part
			obconversion.ReadFile(&templatemol, t_file);
			if (!similarity(targetmol, templatemol, perc))
			{
				if (tmpmols.size() == 0)
				{
					tmpmols.push_back(templatemol);
					std::cout << "templatemol: " << templatemol.GetTitle(false) << std::endl;
					bfneeds[i] -= 1;
				}
				else
				{
					bool find = false;
					for (auto &tmp : tmpmols)
						if (strcmp(tmp.GetTitle(false), templatemol.GetTitle(false)) == 0)
						{ // in fact, strcmp should compaire new.title with bf_names before.
							find = true;
							bfneeds[i] -= 1;
							break;
						}
					if (!find)
					{
						tmpmols.push_back(templatemol);
						std::cout << "templatemol: " << templatemol.GetTitle(false) << std::endl;
						bfneeds[i] -= 1;
					}
				}
			} // !similarity
		} // every bf
		tmgroup.push_back(tmpmols);
	}
	// for every tmpmols
	return tmgroup;
}

// read all templates in tmplt_name
std::map<int, std::vector<OBMol>> buildpocket::readalltmplts(std::set<std::string> tmplt_names,
		OBConversion obconversion, OBMol targetmol, double perc)
{
	const std::map<std::string, BasicFragment> &map =
			BasicFragment::getmap(); // to count bf_num.
	std::string t_file = "template.sdf";
	std::map<int, std::vector<OBMol>> tmpmols;
	std::vector<std::string> tar_out; // target_bf
	for (auto &m : map)
	{
		std::string tar_o = bfdecoupling(&(m.second), &targetmol);
		tar_out.push_back(tar_o);
	}
	int t = 0; // dealing with t_th template.sdf
// read all-sdf.sdf
	std::string ass = getenvpath("DEPACT_DATAPATH")+"all-sdf.sdf";
	std::ifstream ifs(ass.c_str());
	if (!ifs.good())
	{
		std::cout << "sdf failure" << std::endl;
		exit(1);
	}
	bool find = false;
	bool start = true;
	std::ofstream ofs;
	while (true)
	{
		std::string line;
		line.clear();
		getline(ifs, line);
		if (!ifs.good()) break;
		if (line.size() == 0) continue;
		if (tmplt_names.count(line) > 0)
			find = true;
		if (find)
		{
			if (start)
				ofs.open(t_file.c_str());
			else
				ofs.open(t_file.c_str(), std::ios::app);
			ofs << line << std::endl;
			ofs.close();
			if (line != "$$$$")
			{
				if (start)
					start = false;
			}
			else
			{
				start = true;
				find = false;
				OBConversion obconversion;
				OBMol templatemol;
				obconversion.SetInFormat("sdf");
				obconversion.ReadFile(&templatemol, t_file);
				if (!similarity(targetmol, templatemol, perc))
				{
					int bfnum = 0; // count bf_num.
					int i = 0;
					for (auto &m : map)
					{
						std::string o1 = bfdecoupling(&(m.second), &templatemol);
						if (o1 == "Yes" && tar_out[i] == "Yes") bfnum--; // key = -bf_num
						i++;
					} // share the same bf.
					if (tmpmols.count(bfnum) == 0)
					{
						std::vector<OBMol> tpm {templatemol};
						tmpmols.insert(std::make_pair(bfnum, tpm));
					}
					else
						tmpmols[bfnum].push_back(templatemol);
				}
				t++;
				if (t % 100000 == 0)
					std::cout << "dealing with " << t << " templates" << std::endl;
			} // find all common bfs.
		} // find a template.sdf
	} // reading all-sdf.sdf

	return tmpmols;
}

void buildpocket::topsimitmplts(std::vector<OBMol> &tmplt_mols, OBMol targetmol,
		std::vector<std::vector<std::string>> bfns, int needs, double perc)
{
	vector<map<double, vector<string>>> ranks_names(bfns.size()); // bfn<-similar, vector<names>>
	vector<set<string>> bfs_names;
	for (auto bfn : bfns)
	{
		set<string> bfs_n;
		bfs_n.insert(bfn.begin(), bfn.end());
		bfs_names.push_back(bfs_n);
	}
// read all-sdf.sdf
	std::string t_file = "template.sdf";
	std::string ass = getenvpath("DEPACT_DATAPATH")+"all-sdf.sdf";
	std::ifstream ifs(ass.c_str());
	if (!ifs.good())
	{
		std::cout << "sdf failure" << std::endl;
		exit(1);
	}
	bool find = false;
	bool start = true;
	std::ofstream ofs;
	set<int> who; // which bf_group has this tmplt.
	while (true)
	{
		std::string line;
		line.clear();
		getline(ifs, line);
		if (!ifs.good()) break;
		if (line.size() == 0) continue;
		for (int i = 0; i < bfs_names.size(); i++)
			if (bfs_names[i].count(line) > 0)
			{
				find = true;
				who.insert(i);
			}
		if (find)
		{
			if (start)
				ofs.open(t_file.c_str());
			else
				ofs.open(t_file.c_str(), std::ios::app);
			ofs << line << std::endl;
			ofs.close();
			if (line != "$$$$")
			{
				if (start)
					start = false;
			}
			else
			{
				start = true;
				find = false;
				OBConversion obconversion;
				OBMol templatemol;
				obconversion.SetInFormat("sdf");
				obconversion.ReadFile(&templatemol, t_file);
				double value = -similar(targetmol, templatemol);
				if (-value <= perc)
				{
					// test coverability
					for (auto it = who.begin(); it != who.end(); it++)
					{
						int w = *it;
						if (ranks_names[w].count(value) > 0)
							ranks_names[w][value].push_back(templatemol.GetTitle());
						else
						{
							vector<string> names {templatemol.GetTitle()};
							ranks_names[w].insert(make_pair(value, names));
						}
					} // record for every bf who has this tmplt.
				}
				who.clear();
			} // calc. similar
		} // find a template.sdf
	} // reading all-sdf.sdf
	ifs.close();

	// pick topperc.-value names.
	set<string> tmplt_names;
	for (int i = 0; i < bfns.size(); i++)
	{
		auto bfnames = bfns[i];
		auto rank_names = ranks_names[i];
		int pick = needs;
		for (auto it = rank_names.begin(); it != rank_names.end(); it++)
		{
			if (pick >= it->second.size())
			{
				pick -= it->second.size();
				tmplt_names.insert(it->second.begin(), it->second.end());
			}
			else if (pick > 0)
			{
				// randomly choose pick ids in it->second.size()
				for (int id = 0; pick > 0 && id < it->second.size(); id++)
					if (rand() % (it->second.size()-id) < pick)
					{
						tmplt_names.insert(it->second[id]);
						pick--;
					}
				break;
			}
			else break;
		}
		// check: if bf is too large & tmplts are not enough.
		if (pick > 0)
			std::cout << "Warning: bf " << i << " is too large/special to have enough templates." << std::endl;
	}

	// read all-sdf.sdf again to find unredundent tmplt_mols.
	ofstream ofrank("tmplt_similarity.txt");
	std::ifstream ifa(ass.c_str());
	if (!ifa.good())
	{
		std::cout << "sdf failure" << std::endl;
		exit(1);
	}
	find = false;
	start = true;
	while (true)
	{
		std::string line;
		line.clear();
		getline(ifa, line);
		if (!ifa.good()) break;
		if (line.size() == 0) continue;
		if (tmplt_names.count(line) > 0)
			find = true;
		if (find)
		{
			if (start)
				ofs.open(t_file.c_str());
			else
				ofs.open(t_file.c_str(), std::ios::app);
			ofs << line << std::endl;
			ofs.close();
			if (line != "$$$$")
			{
				if (start)
					start = false;
			}
			else
			{
				start = true;
				find = false;
				OBConversion obconversion;
				OBMol templatemol;
				obconversion.SetInFormat("sdf");
				obconversion.ReadFile(&templatemol, t_file);
				tmplt_mols.push_back(templatemol);
				double simi = similar(targetmol, templatemol);
				ofrank << templatemol.GetTitle() << " " << simi << std::endl;
			} // calc. similar
		} // find a template.sdf
	} // reading all-sdf.sdf
	ifa.close();
	ofrank.close();
}

// read tmplt_OBMol from name, and rank them by atm_cover_num. return <-atm_contact_num, vector<id>>
vector<vector<int>> buildpocket::readandrank(std::vector<OBMol> &tmplt_obms,
		std::vector<std::shared_ptr<SubstrAlignment>> &tmplt_ssas,
		std::vector<string> &tmplt_acode, std::set<std::string> tmplt_names,
		OBConversion ob, OBMol targetmol, double perc, int needs, double th)
{
	const std::map<std::string, BasicFragment> &map =
			BasicFragment::getmap(); // to count bf_num.
	std::map<std::string, std::vector<std::shared_ptr<FMMatches>>> fragtargetmatches;
	for (auto &m : map)
		fragtargetmatches[m.first] = findfmmatches(&(m.second), &targetmol, false);
	std::string t_file = "template.sdf";
	std::map<int, set<string>> coverage_code; // -atom_coverage, <code>
	std::map<string, vector<shared_ptr<SubstrAlignment>>> code_obm; // code, <ptr<SA>>
	vector<int> atm_needs(targetmol.NumAtoms(), needs);
// read all-sdf.sdf
	int t = 0; // dealing with t_th template.sdf
	std::string ass = getenvpath("DEPACT_DATAPATH")+"all-sdf.sdf";
	std::ifstream ifs(ass.c_str());
	if (!ifs.good())
	{
		std::cout << "sdf failure" << std::endl;
		exit(1);
	}
	bool find = false;
	bool start = true;
	std::ofstream ofs;
	while (true)
	{
		std::string line;
		line.clear();
		getline(ifs, line);
		if (!ifs.good()) break;
		if (line.size() == 0) continue;
		if (tmplt_names.count(line) > 0)
			find = true;
		if (find)
		{
			if (start)
				ofs.open(t_file.c_str());
			else
				ofs.open(t_file.c_str(), std::ios::app);
			ofs << line << std::endl;
			ofs.close();
			if (line != "$$$$")
			{
				if (start)
					start = false;
			}
			else
			{
				start = true;
				find = false;
				OBConversion obconversion;
				OBMol templatemol;
				obconversion.SetInFormat("sdf");
				obconversion.ReadFile(&templatemol, t_file);
//
//				// test coverability
//				double value = similar(targetmol, templatemol);
//				std::cout << value << std::endl;
//				exit(1);

				if (!similarity(targetmol, templatemol, perc))
				{
					std::map<std::string, std::vector<std::shared_ptr<FMMatches>>> fragtmpltmatches;
					for (auto &m : map)
						fragtmpltmatches[m.first] = findfmmatches(&(m.second), &templatemol, true);
					std::vector<std::shared_ptr<SubstrAlignment>> local_ssas;
					local_ssas = findssas(fragtargetmatches, fragtmpltmatches); // include the combined clique.
					for (auto ssa : local_ssas)
					{
						// calc. atom_coverage(s) = n
						int n = ssa->size();
						// calc. atom_code(s)
						string atm_code(atm_needs.size(), '0');
						for (int i = 0; i < n; ++i)
						{
							auto pa = ssa->alignedpair(i);
							atm_code[pa.first-1] = '1';
							// atom index in openbabel starts from 1
						}
						// renew coverage_code;
						if (coverage_code.count(-n) > 0)
							coverage_code[-n].insert(atm_code);
						else
						{
							set<string> cd {atm_code};
							coverage_code.insert(std::make_pair(-n, cd));
						}
						// renew code_obm;
						if (code_obm.count(atm_code) > 0)
							code_obm[atm_code].push_back(ssa);
						else
						{
							std::vector<std::shared_ptr<SubstrAlignment>> tms {ssa};
							code_obm.insert(std::make_pair(atm_code, tms));
						}
					}
				}
				t++;
				if (t % 100000 == 0)
					std::cout << "dealing with " << t << " templates" << std::endl;
			} // find all common bfs.
		} // find a template.sdf
	} // reading all-sdf.sdf

	int finish = atm_needs.size();
	std::vector<vector<int>> groupids;
	int id = 0;
	for (auto it1 = coverage_code.begin(); it1 != coverage_code.end(); it1++)
	{
		if (finish == 0) break;
		std::vector<int> gp;
		for (auto code : it1->second)
		{
			if (finish == 0) break;
			bool fetch = false;
			for (int i = 0; i < code.size(); i++)
			{
				auto c = code[i];
				if (c == '1' && atm_needs[i] > 0)
				{
					fetch = true;
					break;
				}
			}
			if (fetch)
			{
				int ob_num = ceil(code_obm[code].size() * th);
				for (int i = 0; i < code.size(); i++)
				{
					if (code[i] == '1' && atm_needs[i] > 0)
					{
						atm_needs[i] -= ob_num;
						if (atm_needs[i] <= 0) finish--;
					}
				}
				// fetch := groupids.renew.
				tmplt_acode.push_back(code);
				for (auto &co : code_obm[code])
					tmplt_ssas.push_back(co);
				gp.push_back(id);
				id++;
			}
		} // for all obs(picked) in one code.
		groupids.push_back(gp);
	} // for every code in one atm_coverage.

	return groupids;
}
// fetch the right ssa matches the code
std::vector<std::shared_ptr<SubstrAlignment>> buildpocket::ssaincode(string code, std::vector<std::shared_ptr<SubstrAlignment>> ssas)
{

}
// pick the top th local_ssites and push_back to ssites.
void buildpocket::pickssites(std::vector<TmpltSSAs> &ssas, TargetStruct &tgt, double weight,
		std::vector<std::shared_ptr<Subsite>> &ssites, int pick)
{
	std::vector<std::shared_ptr<Subsite>> local_ss;
	std::map<double, vector<int>> score_rank; // maybe many ss share the same score.
	int id = 0;
	for (auto &ssa : ssas)
	{
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
					local_ss.push_back(ss);
					double tscore = ss->totalscore(weight);
					if (score_rank.count(tscore) > 0)
						score_rank[tscore].push_back(id);
					else
					{
						vector<int> ids {id};
						score_rank.insert(std::make_pair(tscore, ids));
					}
					id++;
				}
			}
			else
				std::cout << "wrong resnum " << residuenumber << " in " << chainid << std::endl;
		}
	} // fetch all local_ss & scores

	for (auto it = score_rank.begin(); it != score_rank.end(); it++)
	{
		if (pick == 0) break;
		for (auto id : it->second)
		{
			ssites.push_back(local_ss[id]); // may bring fault.
			pick--;
			if (pick == 0) break;
		}
	}
}

void buildpocket::makessites(TmpltSSAs ssa, TargetStruct tgt, double weight,
		std::vector<std::shared_ptr<Subsite>> &pre_ssites, int &idx,
		map<int, vector<string>> &covernum_codes,
		map<string, map<double, vector<int>>> &code_ids)
{
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
				string code(tgt.conformer.atomlist.size(), '0'); // empty_atm_code
				for (auto ap : aln.alignedpairs)
				{
					code[ap.first] = '1'; // ap.first has already -1 during ssas->tmpltssas
					codenum--;
				}
				if (covernum_codes.count(codenum) > 0)
					covernum_codes[codenum].push_back(code);
				else
				{
					vector<string> codes {code};
					covernum_codes.insert(make_pair(codenum, codes));
				}
				double tscore = ss->totalscore(weight);
				if (code_ids.count(code) > 0)
				{
					if (code_ids[code].count(tscore) > 0)
					{
						code_ids[code][tscore].push_back(idx);
					}
					else
					{
						vector<int> idxs {idx};
						code_ids[code].insert(make_pair(tscore, idxs));
					}
				}
				else
				{
					map<double, vector<int>> score_idxs;
					vector<int> idxs {idx};
					score_idxs.insert(make_pair(tscore, idxs));
					code_ids.insert(make_pair(code, score_idxs));
				}
				pre_ssites.push_back(ss);
				idx++;
			}
		}
		else
			std::cout << "wrong resnum " << residuenumber << " in " << chainid << std::endl;
	}
}

// pre -> ssites: top_th_ssites
void buildpocket::topssites(std::vector<std::shared_ptr<Subsite>> pre_ssites,
		std::vector<std::shared_ptr<Subsite>> &ssites,
		map<string, map<double, vector<int>>> pcode_ids,
		map<string, map<double, vector<int>>> &code_ids, double th)
{
	// delete tail_code_ids by th.
	for (auto pci = pcode_ids.begin(); pci != pcode_ids.end(); pci++)
	{
		int total_size = 0;
		for (auto sc = pci->second.begin(); sc != pci->second.end(); sc++)
			total_size += sc->second.size();
		int pick = ceil(total_size * th);
		auto d_sc = pci->second.begin();
		for (auto sc = pci->second.begin(); sc != pci->second.end(); sc++)
		{
			pick -= sc->second.size();
			if (pick <= 0)
			{
				d_sc = sc++;
				break;
			}
		}
		pci->second.erase(d_sc, pci->second.end());
	}
	// copy rest_top pre_ssites to ssites. & pcode_ids to code_ids
	int idc = 0;
	for (auto pci = pcode_ids.begin(); pci != pcode_ids.end(); pci++)
	{
		std::map<double, vector<int>> cids;
		for (auto sc = pci->second.begin(); sc != pci->second.end(); sc++)
		{
			vector<int> cis;
			for (auto id : sc->second)
			{
				ssites.push_back(pre_ssites[id]);
				cis.push_back(idc++);
			}
			cids.insert(std::make_pair(sc->first, cis));
		}
		code_ids.insert(std::make_pair(pci->first, cids));
	}
}

// pre_ssites -> unredundent_ssites: Atom_different, otherwise RMSD > 0.3
//void buildpocket::unredundent(std::vector<std::shared_ptr<Subsite>> pre_ssites,
//		std::vector<std::shared_ptr<Subsite>> &ussites,
//		map<string, map<double, vector<int>>> pcode_ids,
//		map<string, map<double, vector<int>>> &ucode_ids, double rmsd)
//{}

bool buildpocket::compactsc(std::vector<NSPgeometry::XYZ> sccrds, Subsite sub_new)
{
	bool compact = true;
	auto krss = sub_new.keyresidues();
	for (auto &krs : krss)
		for (auto &kr : krs)
			for (auto &krcrd : kr.conformer.globalcrd)
			{
				for (auto &scc : sccrds)
				{
					if ((scc-krcrd.second).squarednorm() < DIST2CLASH)
					{
						compact = false;
						goto L1;
					}
				}
			}
	L1: ;
	return compact;
}

void buildpocket::rankformssites(std::vector<std::shared_ptr<Subsite>> &ssites,
		double thread, std::map<double, std::vector<int>> test_scores,
		std::vector<std::shared_ptr<Subsite>> test_ssites)
{
	std::vector<int> best_idx;
	std::cout << "test_ssites.size is: " << test_ssites.size() << std::endl;
	int size = floor(test_ssites.size() * thread);
	for (auto iter = test_scores.begin(); iter != test_scores.end(); iter++)
	{
		if (best_idx.size() > size)
		{
			std::cout << "final ssite's score is " << iter->first << std::endl;
			break;
		}
		if (iter->first >= 0)
		{
			std::cout << "final ssite's score > 0" << std::endl;
			break;
		}
		for (auto i : iter->second)
		{
			auto iter_b = find(best_idx.begin(), best_idx.end(), i);
			if (iter_b == best_idx.end())
				best_idx.push_back(i);
		}
	}
/* test
	std::cout << "test_scores: " << std::endl;
	for (auto i : test_scores)
	{
		std::cout << i.first << " ";
		for (auto s : i.second)
			std::cout << s << " ";
		std::cout << std::endl;
	}
	std::cout << "best top: " << std::endl;
	for (auto i : best_idx)
		std::cout << i << " ";
	std::cout << std::endl;
*/
	for (auto i : best_idx)
		ssites.push_back(test_ssites[i]);

	std::cout << "ssites.size is: " << ssites.size() << std::endl;
	if (ssites.size() == 0)
	{
		std::cout << "there is no ssites" << std::endl;
		exit(1);
	}
}

vector<string> buildpocket::maxscorecodes(std::map<std::string, std::map<double, std::vector<int>>> fcode_ids, string tot_code)
{
	vector<string> max_score_codes;
	int max_score = 0;
	for (auto it = fcode_ids.begin(); it != fcode_ids.end(); it++)
	{
		int score = 0;
		for (int i = 0; i < it->first.size(); i++)
			if (tot_code[i] == '1' && it->first[i] == '1') score++;
		if (score > max_score)
		{
			max_score = score;
			max_score_codes.clear();
			max_score_codes.push_back(it->first);
		}
		else if (score == max_score)
			max_score_codes.push_back(it->first);
	}
	return max_score_codes;
}

int buildpocket::unirandid(std::map<double, std::vector<int>> fci)
{
	int tot_n = 0;
	for (auto it = fci.begin(); it != fci.end(); it++)
		tot_n += it->second.size();
	int r = rand() % tot_n + 1;
	int id = -1;
	for (auto it = fci.begin(); it != fci.end(); it++)
	{
		if (it->second.size() < r)
			r -= it->second.size();
		else
			id = it->second[r-1];
	}
	return id;
}

