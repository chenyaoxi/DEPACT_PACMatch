/*
 * statis_metal.cpp
 *
 *  Created on: 2019年12月9日
 *      Author: yxchen
 */

#include "statis_metal.h"
using namespace myobcode;
using namespace subsitedesign;
using namespace OpenBabel;
using namespace NSPproteinrep;
using namespace statismetal;
#define METALCONTACT 3.5
#define HOHCONTACT 4.0

// judge if is residue
bool statismetal::isresidue(std::string resname)
{
	std::vector<std::string> residues {"GLY", "ALA", "VAL", "LEU",
	"ILE", "MET", "TRP", "PHE", "PRO", "SER", "THR", "CYS", "TYR",
	"ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS"};
	bool isr = false;
	if (std::find(residues.begin(), residues.end(), resname) != residues.end())
		isr = true;
	return isr;
}

std::vector<AAConformersInModel> statismetal::getconfsfromtitlefile(std::string filename, int num_th)
{
	std::vector<std::vector<std::string>> sdfs;
	std::ifstream ifs;
	ifs.open(filename.c_str());
	if (!ifs.good())
	{
		std::cout << "no sdfs_title" << std::endl;
		exit(1);
	}
	while (true)
	{
		std::string line;
		getline(ifs, line);
		if (!ifs.good()) break;
		if (line.size() == 0) continue;
		std::vector<std::string> res;
		std::vector<int> delimposi;
		delimposi.push_back(-1);
		for(int i=1;i<line.size();++i){
			if(line[i]=='_') delimposi.push_back(i);
		}
		for(int i=0;i<delimposi.size()-1;++i){
			res.push_back(line.substr(delimposi[i]+1,delimposi[i+1]-(delimposi[i]+1)));
		}
		sdfs.push_back(res);
	}
	ifs.close();

	std::vector<AAConformersInModel> cims;
	// top to 5000 structures
	int num = 0;
	for (auto &sdf : sdfs)
	{
		const std::string pdbid = sdf[0];
		int modelno = std::stoi(sdf[2]);
		auto reader=readpdbmodel(pdbid, 1, modelno);
		if (!reader) continue; // if pdb is bad: reader == nullptr
		AAConformersInModel cim;
		cim.getconformers(*reader);
		cims.push_back(cim);
		num++;
		if (num >= num_th) break;
	}
	return cims;
}

AAConformersInModel statismetal::getconfsfromtitlefile_ssdf(std::vector<std::string> sdf) // only read single sdf
{
	// top to 5000 structures
	int num = 0;
	const std::string pdbid = sdf[0];
	int modelno = std::stoi(sdf[2]);
	auto reader=readpdbmodel(pdbid, 1, modelno);
	AAConformersInModel cim;
	if (reader) // if pdb is bad: reader == nullptr
		cim.getconformers(*reader);
	return cim;
}

std::vector<AAConformersInModel> statismetal::getconfsfromtitlefile_pdb(std::string filename, int num_th)
{
	std::vector<std::string> pdbs;
	std::ifstream ifs;
	ifs.open(filename.c_str());
	if (!ifs.good())
	{
		std::cout << "no pdbs_title" << std::endl;
		exit(1);
	}
	while (true)
	{
		std::string line;
		getline(ifs, line);
		if (!ifs.good()) break;
		if (line.size() == 0) continue;
		pdbs.push_back(line);
	}
	ifs.close();

	std::vector<AAConformersInModel> cims;
	// top to 5000 structures
	int num = 0;
	for (auto &pdb : pdbs)
	{
		std::cout << pdb << std::endl;
		PdbReader pr;
		pr.readpdb(pdb);
		AAConformersInModel cim;
		cim.getconformers(pr);
		cims.push_back(cim);
		num++;
		if (num >= num_th) break;
	}
	return cims;
}

AAConformersInModel statismetal::getconfsfromtitlefile_spdb(std::string filename)
{
	AAConformersInModel cim;
	int num = 0;
	std::cout << filename << std::endl;
	PdbReader pr;
	pr.readpdb(filename);
	cim.getconformers(pr);
	return cim;
}

std::vector<MetalConf> statismetal::aaconf2metalconf(std::vector<AAConformersInModel> cims)
{
	std::vector<MetalConf> mcs;
	int wr = 0; // water-residue contacts count. to make it not too many.
	for (auto &cim : cims)
	{
		for (auto &cim_chain : cim.conformers)
			for (auto &cim_res : cim_chain)
			{
				if (!cim_res.isnaturalaa() && cim_res.atomlist.size() == 1)
				{
					MetalConf mc;
					auto &crd = cim_res.getglobalcrd();
					mc.readmetal(cim_res.residuename, crd.begin()->second.x_,
							crd.begin()->second.y_, crd.begin()->second.z_);
					if (cim_res.residuename == "HOH")
					{
						for (auto &chain : cim.conformers)
							for (auto &res : chain)
							{
								if (res.residuename == "HOH") continue; // only count water-nonwater contacts
								auto &crds = res.getglobalcrd();
								for (auto iter = crds.begin(); iter != crds.end(); iter++)
								{
									if (iter->first[0] == 'H') continue;
									double dd = iter->second.distance(mc.getmetal().crd);
									if (dd > 0 && dd <= HOHCONTACT)
									{
										if (res.residuename == "CA" || res.residuename == "MG" ||
												res.residuename == "MN" || res.residuename == "ZN")
											mc.addcontact(iter->first, iter->second.x_, iter->second.y_, iter->second.z_);
										else
										{
											string nm;
											nm.push_back(iter->first[0]);
											mc.addcontact(nm, iter->second.x_, iter->second.y_, iter->second.z_);
										}
									}
								}
							}
					}
					else // for metal
					{
						for (auto &chain : cim.conformers)
							for (auto &res : chain)
							{
								auto &crds = res.getglobalcrd();
								for (auto iter = crds.begin(); iter != crds.end(); iter++)
								{
									if (iter->first[0] == 'H') continue;
									double dd = iter->second.distance(mc.getmetal().crd);
									if (dd > 0 && dd <= METALCONTACT)
									{
										if (res.residuename == "CA" || res.residuename == "MG" ||
												res.residuename == "MN" || res.residuename == "ZN")
											mc.addcontact(iter->first, iter->second.x_, iter->second.y_, iter->second.z_);
										else
										{
											string nm;
											nm.push_back(iter->first[0]);
											mc.addcontact(nm, iter->second.x_, iter->second.y_, iter->second.z_);
										}
									}
								}
							}
					}
					mcs.push_back(mc);
				} // for metal & water
			} // every residue
	}
	return mcs;
}

std::vector<MetalConf> statismetal::aaconf2metalconf_single(AAConformersInModel cim)
{
	std::vector<MetalConf> mcs;
	for (auto &cim_chain : cim.conformers)
		for (auto &cim_res : cim_chain)
		{
			if (!cim_res.isnaturalaa() && cim_res.atomlist.size() == 1)
			{
				MetalConf mc;
				auto &crd = cim_res.getglobalcrd();
				mc.readmetal(cim_res.residuename, crd.begin()->second.x_,
						crd.begin()->second.y_, crd.begin()->second.z_);
				if (cim_res.residuename == "HOH")
				{
					for (auto &chain : cim.conformers)
						for (auto &res : chain)
						{
							if (res.residuename == "HOH") continue; // only count water-nonwater contacts
							auto &crds = res.getglobalcrd();
							for (auto iter = crds.begin(); iter != crds.end(); iter++)
							{
								if (iter->first[0] == 'H') continue;
								double dd = iter->second.distance(mc.getmetal().crd);
								if (dd > 0 && dd <= HOHCONTACT)
								{
									if (res.residuename == "CA" || res.residuename == "MG" ||
											res.residuename == "MN" || res.residuename == "ZN")
										mc.addcontact(iter->first, iter->second.x_, iter->second.y_, iter->second.z_);
									else
									{
										string nm;
										nm.push_back(iter->first[0]);
										mc.addcontact(nm, iter->second.x_, iter->second.y_, iter->second.z_);
									}
								}
							}
						}
				}
				else // for metal
				{
					for (auto &chain : cim.conformers)
						for (auto &res : chain)
						{
							auto &crds = res.getglobalcrd();
							for (auto iter = crds.begin(); iter != crds.end(); iter++)
							{
								if (iter->first[0] == 'H') continue;
								double dd = iter->second.distance(mc.getmetal().crd);
								if (dd > 0 && dd <= METALCONTACT)
								{
									if (res.residuename == "CA" || res.residuename == "MG" ||
											res.residuename == "MN" || res.residuename == "ZN")
										mc.addcontact(iter->first, iter->second.x_, iter->second.y_, iter->second.z_);
									else
									{
										string nm;
										nm.push_back(iter->first[0]);
										mc.addcontact(nm, iter->second.x_, iter->second.y_, iter->second.z_);
									}
								}
							}
						}
				}
				mcs.push_back(mc);
			} // for metal & water
		} // every residue


	return mcs;
}

std::vector<MetalConf> statismetal::aaconf2metalconf_spec(std::vector<AAConformersInModel> cims, string metal)
{
	std::vector<MetalConf> mcs;
	int wr = 0; // water-residue contacts count. to make it not too many.
	for (auto &cim : cims)
	{
		for (auto &cim_chain : cim.conformers)
			for (auto &cim_res : cim_chain)
			{
				if (!cim_res.isnaturalaa() && cim_res.atomlist.size() == 1 && cim_res.residuename == metal)
				{
					MetalConf mc;
					auto &crd = cim_res.getglobalcrd();
					mc.readmetal(cim_res.residuename, crd.begin()->second.x_,
							crd.begin()->second.y_, crd.begin()->second.z_);
					if (cim_res.residuename == "HOH")
					{
						for (auto &chain : cim.conformers)
							for (auto &res : chain)
							{
								if (res.residuename == "HOH") continue; // only count water-nonwater contacts
								auto &crds = res.getglobalcrd();
								for (auto iter = crds.begin(); iter != crds.end(); iter++)
								{
									if (iter->first[0] == 'H') continue;
									double bb = iter->second.distance(mc.getmetal().crd);
									if (bb > 0 && bb < HOHCONTACT)
									{
										if (res.residuename == "CA" || res.residuename == "MG" ||
												res.residuename == "MN" || res.residuename == "ZN")
											mc.addcontact(iter->first, iter->second.x_, iter->second.y_, iter->second.z_);
										else
										{
											string nm;
											nm.push_back(iter->first[0]);
											mc.addcontact(nm, iter->second.x_, iter->second.y_, iter->second.z_);
										}
									}
								}
							}
					}
					else
					{
						for (auto &chain : cim.conformers)
							for (auto &res : chain)
							{
								auto &crds = res.getglobalcrd();
								for (auto iter = crds.begin(); iter != crds.end(); iter++)
								{
									if (iter->first[0] == 'H') continue;
									double bb = iter->second.distance(mc.getmetal().crd);
									if (bb > 0 && bb < METALCONTACT)
									{
										if (res.residuename == "CA" || res.residuename == "MG" ||
												res.residuename == "MN" || res.residuename == "ZN")
											mc.addcontact(iter->first, iter->second.x_, iter->second.y_, iter->second.z_);
										else
										{
											string nm;
											nm.push_back(iter->first[0]);
											mc.addcontact(nm, iter->second.x_, iter->second.y_, iter->second.z_);
										}
									}
								}
							}
					}
					mcs.push_back(mc);
				}
			}
	}
	return mcs;
}

std::vector<MetalConf> statismetal::aaconf2metalconf_spec_single(AAConformersInModel cim, string metal)
{
	std::vector<MetalConf> mcs;
	for (auto &cim_chain : cim.conformers)
		for (auto &cim_res : cim_chain)
		{
			if (!cim_res.isnaturalaa() && cim_res.atomlist.size() == 1 && cim_res.residuename == metal)
			{
				MetalConf mc;
				auto &crd = cim_res.getglobalcrd();
				mc.readmetal(cim_res.residuename, crd.begin()->second.x_,
						crd.begin()->second.y_, crd.begin()->second.z_);
				if (cim_res.residuename == "HOH")
				{
					for (auto &chain : cim.conformers)
						for (auto &res : chain)
						{
							if (res.residuename == "HOH") continue; // only count water-nonwater contacts
							auto &crds = res.getglobalcrd();
							for (auto iter = crds.begin(); iter != crds.end(); iter++)
							{
								if (iter->first[0] == 'H') continue;
								double bb = iter->second.distance(mc.getmetal().crd);
								if (bb > 0 && bb < HOHCONTACT)
								{
									if (res.residuename == "CA" || res.residuename == "MG" ||
											res.residuename == "MN" || res.residuename == "ZN")
										mc.addcontact(iter->first, iter->second.x_, iter->second.y_, iter->second.z_);
									else
									{
										string nm;
										nm.push_back(iter->first[0]);
										mc.addcontact(nm, iter->second.x_, iter->second.y_, iter->second.z_);
									}
								}
							}
						}
				}
				else
				{
					for (auto &chain : cim.conformers)
						for (auto &res : chain)
						{
							auto &crds = res.getglobalcrd();
							for (auto iter = crds.begin(); iter != crds.end(); iter++)
							{
								if (iter->first[0] == 'H') continue;
								double bb = iter->second.distance(mc.getmetal().crd);
								if (bb > 0 && bb < METALCONTACT)
								{
									if (res.residuename == "CA" || res.residuename == "MG" ||
											res.residuename == "MN" || res.residuename == "ZN")
										mc.addcontact(iter->first, iter->second.x_, iter->second.y_, iter->second.z_);
									else
									{
										string nm;
										nm.push_back(iter->first[0]);
										mc.addcontact(nm, iter->second.x_, iter->second.y_, iter->second.z_);
									}
								}
							}
						}
				}
				mcs.push_back(mc);
			}
		}
	return mcs;
}

std::map<std::string, ConfStatis> statismetal::metalconf2statis(std::vector<MetalConf> mcs)
{
	std::map<std::string, ConfStatis> maps;
	for (auto &mc : mcs)
	{
		auto metal = mc.getmetal();
		auto cis = mc.getci();
		if (maps.count(metal.atmtype) == 0)
		{
			ConfStatis cs;
			maps.insert(std::make_pair(metal.atmtype, cs));
		}
		// bond statis
		if (metal.atmtype == "HOH")
		{
			for (auto &ci : cis)
				if (ci.atmtype == "O" || ci.atmtype == "N" || ci.atmtype == "S" || ci.atmtype == "C")
					maps[metal.atmtype].bonds.push_back(ci.crd.distance(metal.crd));
				else
					maps[metal.atmtype].wmbonds.push_back(ci.crd.distance(metal.crd));
		}
		else
			for (auto &ci : cis)
				maps[metal.atmtype].bonds.push_back(ci.crd.distance(metal.crd));

		// angle statis
		for (int i = 0; i < cis.size()-1; i++)
			for (int j = i+1; j < cis.size(); j++)
			{
				auto &crdi = cis[i].crd;
				auto &crdj = cis[j].crd;
				const XYZ crda = crdi - metal.crd;
				const XYZ crdb = crdj - metal.crd;
				double cosa = dot(crda, crdb)/(crda.length()*crdb.length());
				maps[metal.atmtype].angles.push_back(acos(cosa)*180/M_PI);
			}
	}
	return maps;
}

std::map<std::string, ConfStatis> statismetal::metalconf2statis(std::vector<MetalConf> mcs, std::string spec_metal)
{
	std::map<std::string, ConfStatis> maps;
	int count = 0;
	std::cout << "mc total size: " << mcs.size() << std::endl;
	for (auto &mc : mcs)
	{

		auto metal = mc.getmetal();
		if (metal.atmtype != spec_metal) continue;
		auto cis = mc.getci();
		if (maps.count(metal.atmtype) == 0)
		{
			ConfStatis cs;
			maps.insert(std::make_pair(metal.atmtype, cs));
		}
		// bond statis
		for (auto &ci : cis)
			maps[metal.atmtype].bonds.push_back(ci.crd.distance(metal.crd));

		// angle statis
		for (int i = 0; i < cis.size()-1; i++)
			for (int j = i+1; j < cis.size(); j++)
			{
				auto &crdi = cis[i].crd;
				auto &crdj = cis[j].crd;
				const XYZ crda = crdi - metal.crd;
				const XYZ crdb = crdj - metal.crd;
				double cosa = dot(crda, crdb)/(crda.length()*crdb.length());
				maps[metal.atmtype].angles.push_back(acos(cosa)*180/M_PI);
			}
	}
	return maps;
}

void statismetal::metalconf2statis_print(std::vector<MetalConf> mcs)
{
	for (auto &mc : mcs)
	{
		auto metal = mc.getmetal();
		auto cis = mc.getci();
		string fb = metal.atmtype + "_bond_statis.txt";
		std::ofstream offb(fb.c_str(), std::ios::app);
		std::string fwb = metal.atmtype + "_metalbond_statis.txt";
		std::ofstream offwb(fwb.c_str(), std::ios::app);
		std::string fa = metal.atmtype + "_angle_statis.txt";
		std::ofstream offa(fa.c_str(), std::ios::app);
		std::string fc = metal.atmtype + "_coord_statis.txt";
		std::ofstream offc(fc.c_str(), std::ios::app);

		// threshold for coordinating num counting
		double th_cn = 0.0;
		if (metal.atmtype == "MG")
			th_cn = 2.7;
		else if (metal.atmtype == "ZN")
			th_cn = 2.7;
		else if (metal.atmtype == "MN")
			th_cn = 2.5;
		else if (metal.atmtype == "CA")
			th_cn = 3.0;
		int cn = 0; // coordinating num.

		// bond statis
		if (metal.atmtype == "HOH")
		{
			for (auto &ci : cis)
				if (ci.atmtype == "O" || ci.atmtype == "N" || ci.atmtype == "S" || ci.atmtype == "C")
				{
					double b = ci.crd.distance(metal.crd);
					offb << b << std::endl;
				}
				else
				{
					double b = ci.crd.distance(metal.crd);
					offwb << b << std::endl;
				}
		}
		else
			for (auto &ci : cis)
			{
				double b = ci.crd.distance(metal.crd);
				offb << b << std::endl;
				if (b <= th_cn) cn++;
			}

		// angle statis
		for (int i = 0; i < cis.size()-1; i++)
			for (int j = i+1; j < cis.size(); j++)
			{
				auto &crdi = cis[i].crd;
				auto &crdj = cis[j].crd;
				const XYZ crda = crdi - metal.crd;
				const XYZ crdb = crdj - metal.crd;
				double cosa = dot(crda, crdb)/(crda.length()*crdb.length());
				double a = acos(cosa)*180/M_PI;
				offa << a << std::endl;
			}

		// coordinating num
		offc << to_string(cn) << std::endl;

		offb.close();
		offwb.close();
		offa.close();
		offc.close();
	}
}
