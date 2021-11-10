/*
 * mmmatches.h
 *
 *  Created on: 2018年8月8日
 *      Author: hyliu
 */

#ifndef MMMATCHES_H_
#define MMMATCHES_H_
#include "move3d.h"
#include "basicfragment.h"
#include "fmmatches.h"
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <memory>
namespace OpenBabel {
class vector3;
}
namespace subsitedesign{
class TmpltSSAs;
}
namespace myobcode {

typedef TwoWayMatches<OpenBabel::OBMol, OpenBabel::OBMol> MMMatches;

/**
 * @brief substructure alignment between two ligand molecules
 */
class SubstrAlignment {
public:
	SubstrAlignment(const FMMatches *fm1, const FMMatches *fm2);
	SubstrAlignment(const SubstrAlignment &ssa1, const SubstrAlignment &ssa2);
	MMMatches * mmmatches() const {
		return mmmatches_.get();
	}
	subsitedesign::Move3D *move3d() const {
		return move3d_.get();
	}
	double rmsd() const {
		return rmsd_;
	}
	const OpenBabel::OBMol *mol1() const {
		return mmmatches_->obj1;
	}
	const OpenBabel::OBMol *mol2() const {
		return mmmatches_->obj2;
	}
	int matchof1in2(int idx1) const {
		return mmmatches_->o1matchino2(idx1);
	}
	int matchof2in1(int idx2) const {
		return mmmatches_->o2matchino1(idx2);
	}
	std::pair<int, int> alignedpair(int idx) const {
		return std::make_pair(mmmatches_->ithmatchino1(idx),
				mmmatches_->ithmatchino2(idx));
	}
/*	const std::vector<int> &specialatoms() const {
		return specialatoms_;
	}
	const std::vector<std::string> &sptypes() const {
		return sptypes_;
	}
	bool atomspecial(int idx) const {
		for(auto i:specialatoms_) if(i==idx) return true;
		return  false;
	}
	std::string sptype(int idx) const {
		for(int i=0;i<specialatoms_.size();++i){
			if(specialatoms_[i]==idx) return sptypes_[i];
		}
		return "";
	}*/
	int size() const {
		return mmmatches_->nmatches();
	}
private:
	std::shared_ptr<MMMatches> mmmatches_ { nullptr };
	//std::vector<int> specialatoms_;
	//std::vector<std::string> sptypes_;
	std::shared_ptr<subsitedesign::Move3D> move3d_ { nullptr };
	double rmsd_ { 100000 };
};
/**
 * matches between two ligand molecules based on two matches to the same basic fragments
 */
std::shared_ptr<MMMatches> fm2mmmatches(const FMMatches *fm1,
		const FMMatches *fm2);
/**
 * determine the rigid tansformation and RMSD for a match (alignment)
 */
std::shared_ptr<subsitedesign::Move3D> fitm2to1(MMMatches *mmmatches, double *rmsd = nullptr);
OpenBabel::vector3 center(const std::vector<OpenBabel::vector3> &atoms);
subsitedesign::Move3D fitatoms(const std::vector<OpenBabel::vector3> & atms1,
		const std::vector<OpenBabel::vector3> &atms2, double *rmsd = nullptr);
void moveobmol(OpenBabel::OBMol *mol, const subsitedesign::Move3D *move3d);
/**
 * @brief determines compatibility of two substructure alignment
 */
// if ssa1 & ssa2 are compatible
bool compatible(const SubstrAlignment &ssa1, SubstrAlignment &ssa2);
// compared with compatible, combinable required ssa1 != ssa2.
bool combinable(const SubstrAlignment &ssa1, SubstrAlignment &ssa2);
bool different(SubstrAlignment ssa, std::vector<std::shared_ptr<SubstrAlignment>> ssas);
/**
 * @brief find substructure alignments between a target and a template ligand
 *
 * First alignments based on basic fragments are obtained
 * Then the maximal clique approach is used to find the extended substructure alignments (ssas).
 */
std::vector<std::shared_ptr<SubstrAlignment>> findssas(
		const std::map<std::string,std::vector<std::shared_ptr<FMMatches>>> &fragtargetmatches,
		const std::map<std::string,std::vector<std::shared_ptr<FMMatches>>> &fragtmpltmatches);
/**
 * @brief construct extended ssas (maximal cliques) from the input ssas
 *
 */

// extendssas: try to find max_clique for combining ssas
void extendssas(std::vector<std::shared_ptr<SubstrAlignment>> &ssas);

std::shared_ptr<subsitedesign::TmpltSSAs> maketmpltssas(
		const std::vector<std::shared_ptr<SubstrAlignment>> &ssas,const std::vector<std::string> &tmplatomtypes);
}

#endif /* MMMATCHES_H_ */
