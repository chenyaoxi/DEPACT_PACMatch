/*
 *
 * mmmatches.cpp
 *
 *  Created on: 2018年8月8日
 *      Author: hyliu
 */

#include "mmmatches.h"
#include "openbabel/math/vector3.h"
#include "tmpltssas.h"
#include <set>
#include <cassert>
using namespace OpenBabel;
using namespace myobcode;
using namespace subsitedesign;
static double rmsdmax(int natoms){
	if(natoms<=4) return 0.1;
	else if(natoms <=8) return 0.25;
//	else if(natoms <=16) return 0.4;
	else return 0.5;
}
static vector3 tovector3(const std::vector<double> &x) {
	vector3 res(x[0], x[1], x[2]);
	return res;
}
static std::vector<double> tostdvector(const vector3& x) {
	std::vector<double> res(3);
	for (int m = 0; m < 3; ++m)
		res[m] = x[m];
	return res;
}

bool myobcode::different(SubstrAlignment ssa1, std::vector<std::shared_ptr<SubstrAlignment>> ssas)
{
	for (auto pssa2 : ssas)
	{
		auto ssa2 = *pssa2;
		int n1 = ssa1.size();
		int n2 = ssa2.size();
		if (n1 == n2)
		{
			std::set<std::pair<int, int>> p1s, p2s;
			bool same = true;
			for (int i = 0; i < n1; ++i)
			{
				std::pair<int, int> p1 = ssa1.alignedpair(i);
				std::pair<int, int> p2 = ssa2.alignedpair(i);
				p1s.insert(p1);
				p2s.insert(p2);
			}
			auto it2 = p2s.begin();
			for (auto it1 = p1s.begin(); it1 != p1s.end(); it1++)
			{
				if (it1->first != it2->first || it2->second != it2->second)
				{
					same = false;
					break;
				}
				it2++;
			}
			if (same) return false;
		} // have same size
	} // for every ssa in ssas

	return true;
}

/*
this version is too slow comparing with that in extendssas.cpp
void myobcode::extendssas(std::vector<std::shared_ptr<SubstrAlignment>> &ssas)
{
	std::cout << "this is my extendssas" << std::endl;
	// ver1. save all alignments: cost time and memory
	bool hasnew = true;
	while(hasnew)
	{
		hasnew = false;
		for (int i = 0; i < ssas.size()-1; i++)
			for (int j = i+1; j < ssas.size(); j++)
			{
				if (combinable(*(ssas[i]), *(ssas[j])))
				{
					auto nssa = std::shared_ptr < SubstrAlignment
							> (new SubstrAlignment(*(ssas[i]), *(ssas[j])));
					if (different(*nssa, ssas))
					{
						ssas.push_back(nssa);
						hasnew = true;
					}
				}
			}
	}


/*
	// ver2. only extend the maximun cilques.
	std::map<int, std::set<int>> map_ssas; // id_ssas, <combinable_id_ssas, including itself for convinence>
	for (int i = 0; i < ssas.size(); i++)
	{
		std::set<int> setid;
		setid.insert(i);
		map_ssas.insert(std::make_pair(i, setid));
	}
	for (int i = 0; i < ssas.size()-1; i++)
		for (int j = 1; j < ssas.size(); j++)
			if (combinable(*(ssas[i]), *(ssas[j])))
			{
					map_ssas[i].insert(j);
					map_ssas[j].insert(i);
			}
}
*/

std::vector<std::shared_ptr<SubstrAlignment>> myobcode::findssas(
		const std::map<std::string, std::vector<std::shared_ptr<FMMatches>>>&fragtargetmatches,
const std::map<std::string,std::vector<std::shared_ptr<FMMatches>>> &fragtmpltmatches) {
	const std::map<std::string, BasicFragment> &smmap =
	BasicFragment::getmap();
//	double rmsdmax=0.3;
	std::vector<std::string> smartsid;
	std::vector<std::shared_ptr<SubstrAlignment>> ssas;
	for(auto &m:smmap) {
		const std::vector<std::shared_ptr<FMMatches>> &fmmtarget=fragtargetmatches.at(m.first);
		if(fmmtarget.empty()) continue;
		const std::vector<std::shared_ptr<FMMatches>> &fmmtemplate =
		fragtmpltmatches.at(m.first);
		for (auto mtgt : fmmtarget) {
			for (auto mtmplt : fmmtemplate) {
				auto nssa = std::shared_ptr < SubstrAlignment
				> (new SubstrAlignment(mtgt.get(), mtmplt.get()));
				if (nssa->rmsd() > rmsdmax(nssa->size()))
				continue;
				/*			int idx = 0;
				 bool redundant = false;
				 for (auto &ossa : ssas) {
				 if (equivalent(*nssa, *ossa)) {
				 ossa = std::shared_ptr < SubstrAlignment
				 > (new SubstrAlignment(*nssa, *ossa));
				 smartsid[idx] = smartsid[idx] + "_M_" + m.first;
				 redundant = true;
				 break;
				 }
				 idx++;
				 }
				 if (!redundant) {*/
				ssas.push_back(nssa);
				smartsid.push_back(m.first);
//					idx++;}
			}
		}
	}
	extendssas(ssas); // it's necessary to bring largest ssas.

	return ssas;
}
/*
void myobcode::mergeequivalentssas(
		std::vector<std::shared_ptr<SubstrAlignment>> &ssas) {
	std::vector<std::shared_ptr<SubstrAlignment>> oldssas = ssas;
	ssas.resize(1);
	for (int i = 1; i < oldssas.size(); ++i) {
		auto nssa = oldssas[i];
		bool redundant = false;
		for (auto &ossa : ssas) {
			if (equivalent(*nssa, *ossa)) {
				ossa = std::shared_ptr < SubstrAlignment
						> (new SubstrAlignment(*nssa, *ossa));
				redundant = true;
				break;
			}
		}
		if (!redundant) {
			ssas.push_back(nssa);
		}
	}
}*/
SubstrAlignment::SubstrAlignment(const FMMatches *fm1, const FMMatches *fm2) {
	mmmatches_ = fm2mmmatches(fm1, fm2);
	move3d_ = fitm2to1(mmmatches_.get(), &rmsd_);
//	specialatoms_ = fm1->obj1->spatoms;
//	sptypes_ = fm1->obj1->sptypes;
}
SubstrAlignment::SubstrAlignment(const SubstrAlignment &ssa1,
		const SubstrAlignment &ssa2) {
	mmmatches_ = std::shared_ptr < MMMatches > (new MMMatches());
	*mmmatches_ = *(ssa1.mmmatches());
	const std::vector<int> & ssa2seq1 = ssa2.mmmatches()->seq1;
	const std::vector<int> &ssa2seq2 = ssa2.mmmatches()->seq2;
	//const std::vector<int> ssa2spatms = ssa2.specialatoms();
	std::vector<int> &seq1 = mmmatches_->seq1;
	std::vector<int> &seq2 = mmmatches_->seq2;
	std::map<int, int> &map12 = mmmatches_->map12;
	std::map<int, int> &map21 = mmmatches_->map21;
//	specialatoms_ = ssa1.specialatoms();
//	sptypes_ = ssa1.sptypes();
	for (int i = 0; i < ssa2seq1.size(); ++i) {
		bool found = false;
		for (int j = 0; j < seq1.size(); ++j) {
			if (ssa2seq1.at(i) == seq1[j]) {
				assert(ssa2seq2.at(i) == seq2[j]);
				found = true;
				break;
			}
		}
		if (!found) {
			seq1.push_back(ssa2seq1.at(i));
			seq2.push_back(ssa2seq2.at(i));
			map12[ssa2seq1.at(i)] = ssa2seq2.at(i);
			map21[ssa2seq2.at(i)] = ssa2seq1.at(i);
/*			for (auto s : ssa2spatms)
				if (s == i) {
					specialatoms_.push_back(seq1.size() - 1);
					sptypes_.push_back(ssa2.sptype(s));
					assert(!sptypes_.back().empty());
				}*/
		}
	}
	move3d_ = fitm2to1(mmmatches_.get(), &rmsd_);
}
bool myobcode::compatible(const SubstrAlignment &ssa1, SubstrAlignment &ssa2) {
	bool result = true;
	result = (ssa1.mol1() == ssa2.mol1()) && (ssa1.mol2() == ssa2.mol2());
	if (!result)
		return result;
	int n1 = ssa1.size();
	int n2 = ssa2.size();
	std::set<std::pair<int, int>> mergedpairs;
	for (int i = 0; i < n1; ++i) {
		std::pair<int, int> p1 = ssa1.alignedpair(i);
		mergedpairs.insert(p1);
		for (int j = 0; j < n2; ++j) {
			std::pair<int, int> p2 = ssa2.alignedpair(j);
			if (p1.first == p2.first) {
				if (p1.second != p2.second) {
					result = false;
					break;
				}
			}
			mergedpairs.insert(p2);
		}
		if (!result)
			return result;
	}
	if (mergedpairs.size() == n1)
		return result;
	std::vector<vector3> atms1;
	std::vector<vector3> atms2;
	for (auto &p : mergedpairs) {
		atms1.push_back(ssa1.mol1()->GetAtom(p.first)->GetVector());
		atms2.push_back(ssa1.mol2()->GetAtom(p.second)->GetVector());
	}
	double rmsd;
	fitatoms(atms1, atms2, &rmsd);
	if (rmsd > rmsdmax(atms1.size()))
		result = false;
	return result;
}

bool myobcode::combinable(const SubstrAlignment &ssa1, SubstrAlignment &ssa2) {
	bool result = true;
	result = (ssa1.mol1() == ssa2.mol1()) && (ssa1.mol2() == ssa2.mol2());
	if (!result)
		return result;
	int n1 = ssa1.size();
	int n2 = ssa2.size();
	std::set<std::pair<int, int>> mergedpairs;
	for (int i = 0; i < n1; ++i) {
		std::pair<int, int> p1 = ssa1.alignedpair(i);
		mergedpairs.insert(p1);
		for (int j = 0; j < n2; ++j) {
			std::pair<int, int> p2 = ssa2.alignedpair(j);
			if (p1.first == p2.first) {
				if (p1.second != p2.second) {
					result = false;
					break;
				}
			}
			mergedpairs.insert(p2);
		}
		if (!result)
			return result;
	}
	if (mergedpairs.size() <= n1 || mergedpairs.size() <= n2)
	{
		result = false; // if merge donot become larger: false.
		return result;
	}
	std::vector<vector3> atms1;
	std::vector<vector3> atms2;
	for (auto &p : mergedpairs) {
		atms1.push_back(ssa1.mol1()->GetAtom(p.first)->GetVector());
		atms2.push_back(ssa1.mol2()->GetAtom(p.second)->GetVector());
	}
	double rmsd;
	fitatoms(atms1, atms2, &rmsd);
	if (rmsd > rmsdmax(atms1.size()))
		result = false;
	return result;
}
/*bool myobcode::equivalent(const SubstrAlignment &ssa1, SubstrAlignment &ssa2) {
	bool result = true;
	result =
			(ssa1.mol1() == ssa2.mol1())
					&& (ssa2.mol2() == ssa2.mol2()
							&& (ssa1.specialatoms().size()
									== ssa2.specialatoms().size()));
	if (!result)
		return result;
	std::set<std::pair<int, int>> sp1;
	for (int i = 0; i < ssa1.specialatoms().size(); ++i)
		sp1.insert(ssa1.alignedpair(ssa1.specialatoms().at(i)));
	for (int i = 0; i < ssa2.specialatoms().size(); ++i) {
		if (sp1.find(ssa2.alignedpair(ssa2.specialatoms().at(i)))
				== sp1.end()) {
			return false;
		}
	}
	return compatible(ssa1, ssa2);
}*/
std::shared_ptr<MMMatches> myobcode::fm2mmmatches(const FMMatches *fm1,
		const FMMatches *fm2) {
	std::vector<int> list1, list2;
	for (int i = 0; i < fm1->nmatches(); ++i) {
		list1.push_back(fm1->map12.at(i));
		list2.push_back(fm2->map12.at(i));
	}
	return std::shared_ptr < MMMatches
			> (new MMMatches(fm1->obj2, fm2->obj2, list1, list2));
}
std::shared_ptr<Move3D> myobcode::fitm2to1(MMMatches *mmmatches, double *rmsd) {
	const OBMol *mol1 = mmmatches->obj1;
	const OBMol *mol2 = mmmatches->obj2;
	std::vector<vector3> atms1;
	std::vector<vector3> atms2;
	for (int i = 0; i < mmmatches->nmatches(); ++i) {
		atms1.push_back(mol1->GetAtom(mmmatches->ithmatchino1(i))->GetVector());
		atms2.push_back(mol2->GetAtom(mmmatches->ithmatchino2(i))->GetVector());
	}
	return std::shared_ptr < Move3D > (new Move3D(fitatoms(atms1, atms2, rmsd)));
}

void myobcode::moveobmol(OBMol *mol, const Move3D *move3d) {
	vector3 trans = tovector3(move3d->tvec);
	double rmatrix[3][3];
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j)
			rmatrix[i][j] = move3d->rmatrix[i][j];
	}
	mol->Rotate(rmatrix);
	mol->Translate(trans);
}
vector3 myobcode::center(const std::vector<vector3> &atoms) {
	vector3 c(0.0, 0.0, 0.0);
	for (auto &a : atoms)
		c = c + a;
	c = c * (1.0 / (double) atoms.size());
	return c;
}
#include <openbabel/obutil.h>

Move3D myobcode::fitatoms(const std::vector<OpenBabel::vector3> & atms1,
		const std::vector<OpenBabel::vector3> &atms2, double *rmsd) {
	vector3 c1 = center(atms1);
	vector3 c2 = center(atms2);
	int natoms = atms1.size();
	double *r1 = new double[natoms * 3];
	double *r2 = new double[natoms * 3];
	for (int i = 0; i < natoms; ++i) {
		for (int m = 0; m < 3; ++m) {
			r1[3 * i + m] = atms1[i][m] - c1[m];
			r2[3 * i + m] = atms2[i][m] - c2[m];
		}
	}
	Move3D res;
	std::vector<double> &trans = res.tvec;
	double (*rmatrix)[3] = res.rmatrix;
	qtrfit(r1, r2, (unsigned int) natoms, rmatrix);
	for (int i = 0; i < 3; ++i) {
		trans[i] = c1[i]
				- (rmatrix[i][0] * c2[0] + rmatrix[i][1] * c2[1]
						+ rmatrix[i][2] * c2[2]);
	}
	if (rmsd) {
		std::vector<vector3> atms2m;
		for (auto &a : atms2) {
			atms2m.push_back(tovector3(res.move(tostdvector(a))));
		}
		double d2 = 0.0;
		for (int i = 0; i < natoms; ++i) {
			d2 += (atms1[i] - atms2m[i]).length_2();
		}
		*rmsd = sqrt(d2);
	}
	delete[] r1, r2;
	return res;
}
using namespace subsitedesign;
std::shared_ptr<TmpltSSAs> myobcode::maketmpltssas(
		const std::vector<std::shared_ptr<SubstrAlignment>> &ssas,
		const std::vector<std::string> &tmpltatomtypes) {
	auto res = std::shared_ptr < TmpltSSAs > (new TmpltSSAs());
	res->tmpltname = std::string(ssas[0]->mol2()->GetTitle());
	res->targetname = std::string(ssas[0]->mol1()->GetTitle());
	res->atomtypes=tmpltatomtypes;
	for (auto p : ssas) {
		res->alignments.push_back(TmpltSSAs::Alignment());
		TmpltSSAs::Alignment &algn = res->alignments.back();
		int n = p->size();
		for (int i = 0; i < n; ++i) {
			auto pa = p->alignedpair(i);
			// atom index in openbabel starts from 1
			algn.alignedpairs.insert(
					std::make_pair(pa.first - 1, pa.second - 1));
		}
		algn.rmsd = p->rmsd();
		algn.move3d = *(p->move3d());
	}
	return res;
}
