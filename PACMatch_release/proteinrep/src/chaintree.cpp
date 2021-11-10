/*
 * backbonegeometry.cpp
 *
 *  Created on: 2016年11月17日
 *      Author: hyliu
 */
#include "proteinrep/chaintree.h"
#include "proteinrep/idealgeometries.h"
#include "proteinrep/pdbrecord.h"

using namespace NSPproteinrep;
std::shared_ptr<ChainTreeTopology> NSPproteinrep::makeMainChainHeavyTopo(
		unsigned int length, unsigned int chainid) {
	typedef typename ChainTreeTopology::AtomKey AtomKey;
	std::shared_ptr<ChainTreeTopology> toposp(
			new ChainTreeTopology(length, chainid,
					BackBoneData::atomnames[BackBoneData::N]));
	ChainTreeTopology *topo = toposp.get();
	for (unsigned int i = 0; i < length; ++i) {
		AtomKey n_k = topo->genKey(i, BackBoneData::atomnames[BackBoneData::N]);
		AtomKey ca_k = topo->genKey(i,
				BackBoneData::atomnames[BackBoneData::CA]);
		AtomKey c_k = topo->genKey(i, BackBoneData::atomnames[BackBoneData::C]);
		AtomKey o_k = topo->genKey(i, BackBoneData::atomnames[BackBoneData::O]);
		if (i > 0) {
			AtomKey p_k = topo->genKey(i - 1,
					BackBoneData::atomnames[BackBoneData::C]);
			topo->attachAtom(p_k, n_k, true);
		}
		topo->attachAtom(n_k, ca_k, true);
		topo->attachAtom(ca_k, c_k, true);
		topo->attachAtom(c_k, o_k, false);
	}
	return toposp;
}

void ChainTreeCrd::writePDB(std::ostream &os) const {
	int atomid=1;
	for(auto & c:crdmap_){
		PdbRecord record=make_pdbrecord<AtomKeyTypeA>(c.first,c.second,atomid++,1);
		if(record.residuename == "ANY") record.residuename="GLY";
		os<<record.toString() <<std::endl;
	}

}
void ChainTreeCrd::initWithIdealIntCrd() {
	IdealGeometries & idg = IdealGeometries::getGlobalInstance();
	auto & subtrees = topo()->subtreemap();
	for (auto & st : subtrees) {
		setIdealIntCrd(st.first);
	}
	for (auto & st : subtrees) {
		const AtomKey & la = st.first;
		NSPgeometry::IntCrd & ic = intcrdmap_.at(la);
		if (NSPgeometry::IntCrd::validTorsion(ic.torsion))
			continue;
		setIdealRelativeTorsion(la);
	}
}
void ChainTreeCrd::setBackBoneTorsionAt(unsigned int posi,const BackBoneTorsions & tor){
	setPhi(posi,tor.phi);
	setPsi(posi,tor.psi);
	setOmega(posi,tor.omega);
}
void ChainTreeCrd::setBackBoneTorsions(const std::map<unsigned int,BackBoneTorsions> &tors){
	for(auto & t:tors) {
		setBackBoneTorsionAt(t.first,t.second);
	}
}
void ChainTreeCrd::setIdealIntCrd(const AtomKey &la) {
	IdealGeometries & idg = IdealGeometries::getGlobalInstance();
	auto & aijk = topo()->aijk(la);
	NSPgeometry::IntCrd & ic = intcrdmap_.at(la);
	double d = idg.idealLength(ChainTreeTopology::atomName(la),
			ChainTreeTopology::atomName(aijk.ka));
	if (NSPgeometry::IntCrd::validDistance(d))
		ic.distance = d;
	double a = idg.idealAngle(ChainTreeTopology::atomName(la),
			ChainTreeTopology::atomName(aijk.ka),
			ChainTreeTopology::atomName(aijk.ja));
	if (NSPgeometry::IntCrd::validAngle(a))
		ic.angle = a;
	double t = idg.idealTorsion(ChainTreeTopology::atomName(la),
			ChainTreeTopology::atomName(aijk.ka),
			ChainTreeTopology::atomName(aijk.ja),
			ChainTreeTopology::atomName(aijk.ia));
	if (NSPgeometry::IntCrd::validTorsion(t))
		ic.torsion = t;
}
void ChainTreeCrd::setIdealRelativeTorsion(const AtomKey &la) {
	IdealGeometries & idg = IdealGeometries::getGlobalInstance();
	NSPgeometry::IntCrd & ic = intcrdmap_.at(la);
	auto & aijk = topo()->aijk(la);
	auto kachildren = topotree_->branch(aijk.ka)->children();
	for (auto & b : kachildren) {
		NSPgeometry::IntCrd & icn = intcrdmap_.at(b.first);
		if (NSPgeometry::IntCrd::validTorsion(icn.torsion)) {
			double dt = idg.idealRelativeTorsion(
					ChainTreeTopology::atomName(aijk.ja),
					ChainTreeTopology::atomName(aijk.ka),
					ChainTreeTopology::atomName(b.first),
					ChainTreeTopology::atomName(la));
			if (NSPgeometry::IntCrd::validTorsion(dt)) {
				ic.torsion = icn.torsion + dt;
				break;
			}
		}
	}
}
void ChainTreeCrd::copyCrdFromSites(std::vector<BackBoneSite>::const_iterator iter,
		bool firstcis){
	std::map<AtomKey,NSPgeometry::XYZ> crd;
	for (int i=0; i<length();++i){
		AtomKey n_k=topo()->genKey(i,BackBoneData::atomnames[BackBoneData::N]);
		AtomKey ca_k=topo()->genKey(i,BackBoneData::atomnames[BackBoneData::CA]);
		AtomKey c_k=topo()->genKey(i,BackBoneData::atomnames[BackBoneData::C]);
		AtomKey o_k=topo()->genKey(i,BackBoneData::atomnames[BackBoneData::O]);
		auto bs=iter+i;
		crd.insert(std::make_pair(n_k,bs->getcrd(BackBoneSite::NCRD)));
		crd.insert(std::make_pair(ca_k,bs->getcrd(BackBoneSite::CACRD)));
		crd.insert(std::make_pair(c_k,bs->getcrd(BackBoneSite::CCRD)));
		crd.insert(std::make_pair(o_k,bs->getcrd(BackBoneSite::OCRD)));
	}
	copyCrdMap(crd);
	double deg=3.14159265358979323846 / 180.0;
	setPhi(0,iter->phi()*deg);
	setPsi(0,iter->psi()*deg);
	if(firstcis) setOmega(0,0.0);
	else setOmega(0,3.14159265358979323846);
	setPsi(length()-1,(iter+length()-1)->psi()*deg);
}
void ChainTreeCrd::initCrdPosi0(const std::map<std::string,NSPgeometry::XYZ> & crd0){
	std::map<AtomKey,NSPgeometry::XYZ> crd;
	auto n_it=crd0.find(BackBoneData::atomnames[BackBoneData::N]);
	if( n_it != crd0.end()) {
		for(auto & c:crd0) {
			AtomKey k=topo()->genKey(0,c.first);
			crd.insert(std::make_pair(k,c.second));
		}
	} else {
		IdealGeometries & idg = IdealGeometries::getGlobalInstance();
		AtomKey n_k=topo()->genKey(0,BackBoneData::atomnames[BackBoneData::N]);
		AtomKey ca_k=topo()->genKey(0,BackBoneData::atomnames[BackBoneData::CA]);
		AtomKey c_k=topo()->genKey(0,BackBoneData::atomnames[BackBoneData::C]);
		NSPgeometry::XYZ rn(0.0,0.0,0.0);
		NSPgeometry::XYZ rca=NSPgeometry::InternaltoXYZ(rn,
				idg.idealLength(BackBoneData::atomnames[BackBoneData::N],BackBoneData::atomnames[BackBoneData::CA]));
		NSPgeometry::XYZ rc=NSPgeometry::InternaltoXYZ(rca,rn,
						idg.idealLength(BackBoneData::atomnames[BackBoneData::CA],BackBoneData::atomnames[BackBoneData::C]),
						idg.idealAngle(BackBoneData::atomnames[BackBoneData::N],BackBoneData::atomnames[BackBoneData::CA],
								BackBoneData::atomnames[BackBoneData::C]));
		crd.insert(std::make_pair(n_k,rn));
		crd.insert(std::make_pair(ca_k,rca));
		crd.insert(std::make_pair(c_k,rc));
	}
	copyCrdMap(crd);
}
void ChainTreeCrd::writeIntCrd(std::ostream &os) const {
	auto & subtrees = topo()->subtreemap();
	double rad=180.0/3.14159265;
	for (auto & st : subtrees) {
		const AtomKey & la = st.first;
		auto & ic=intcrdmap_.at(la);
		auto & aijk = topo()->aijk(la);
		std::cout << AtomKeyTypeA::posiNumber(la) <<" " <<
				  AtomKeyTypeA::atomName(la) << " " <<AtomKeyTypeA::atomName(aijk.ka) <<" "
	               << ic.distance <<" "<< AtomKeyTypeA::atomName(aijk.ja) <<" "<< ic.angle*rad
				   <<" " <<AtomKeyTypeA::atomName(aijk.ia) <<" "<<ic.torsion*rad <<std::endl;
	}
}
