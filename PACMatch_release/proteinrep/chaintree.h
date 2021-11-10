/*
 * backbone.h
 *
 *  Created on: 2016年11月17日
 *      Author: hyliu
 */

#ifndef PROTEINREP_CHAINTREE_H_
#define PROTEINREP_CHAINTREE_H_
#include <backbone/backbonedata.h>
#include <backbone/backbonesite.h>
#include "geometry/crdtree.h"
#include "proteinrep/intatomkey.h"
#include <iostream>
#include <memory>
namespace NSPproteinrep {
class ChainTreeTopology : public NSPgeometry::TopoTree<typename AtomKeyTypeA::Key> {
public:
	typedef typename AtomKeyTypeA::Key AtomKey;
	AtomKey genKey(unsigned int posi, const std::string & atomname) const {
		assert(posi >= 0 && posi < length_);
		return AtomKeyTypeA::genKey(posi,atomname,chainid_);
	}

	ChainTreeTopology(unsigned int l, unsigned int chainid,
			const std::string &atomname0):
				TopoTree<typename AtomKeyTypeA::Key>(
						AtomKeyTypeA::genKey(0,atomname0,chainid)),
						length_(l),chainid_(chainid){;}
	unsigned int length() const {return length_;}
	unsigned int chainid() const {return chainid_;}
	static std::string atomName(const AtomKey &key) {
		return AtomKeyTypeA::atomName(key);
	}
	static unsigned int posiNumber(const AtomKey &key) {
		return AtomKeyTypeA::posiNumber(key);
	}
	std::vector<std::pair<unsigned int,std::string>> atomlist() const {
		std::vector<std::pair<unsigned int,std::string>> list;
		const auto & subtrees=subtreemap();
		for(const auto & st:subtrees) {
			list.push_back(std::make_pair(posiNumber(st.first),atomName(st.first)));
		}
		return list;
	}
	template<typename DATA>
	std::map<AtomKey,DATA> remapAtomData(const std::map<std::pair<unsigned int,std::string>,DATA> & data){
		std::map<AtomKey,DATA> res;
		for(auto &d:data) {
			AtomKey k=genKey(d.first.first,d.first.second);
			res.insert(std::make_pair(k,d.second));
		}
		return std::move(res);
	}
private:
	unsigned int length_{0};
	unsigned int chainid_{0};
};
std::shared_ptr<ChainTreeTopology> makeMainChainHeavyTopo(
		unsigned int length,unsigned int chainid=0);

class ChainTreeCrd:  public NSPgeometry::CrdTree<typename ChainTreeTopology::AtomKey>  {
public:
	typedef ChainTreeTopology Topology;
	typedef typename Topology::AtomKey AtomKey;
	typedef NSPgeometry::CrdTree<AtomKey> Coordinate;
	struct BackBoneTorsions {
		BackBoneTorsions() {
			;
		}
		BackBoneTorsions(double Phi, double Psi, double Omega = 3.14159265358979323846) {

			phi = Phi;
			psi = Psi;
			omega = Omega;
		}
		std::string toString() {
			double rad=180.0/3.14159265358979323846;
			return  "(" + std::to_string(shiftangle(phi*rad)) +","
					+ std::to_string(shiftangle(psi*rad)) +"," +std::to_string(shiftangle(omega*rad)) +")";
		}
		double phi { 1000.0 };  //phi,psi initialized to invalid values according to IntCrd::validTorsion
		double psi { 1000.0 };
		double omega { 3.14159265358979323846 };
		double shiftangle(double angle, double min=-180.0, double max=180.0){
			if(angle <min) angle +=360.0;
			if(angle >max) angle -=360.0;
			return angle;
		}
	};
	ChainTreeCrd(ChainTreeTopology * topo):
		Coordinate(static_cast<NSPgeometry::TopoTree<typename AtomKeyTypeA::Key> *>(topo)) {;}
	ChainTreeTopology *topo() {return static_cast<ChainTreeTopology *> (topotree_);}
	const ChainTreeTopology *topo() const {return static_cast<const ChainTreeTopology *> (topotree_);}
	void initWithIdealIntCrd();
	void copyCrdFromSites(std::vector<BackBoneSite>::const_iterator iter, bool firstcis=false);
	void setIdealIntCrd(const AtomKey &key);
	void setIdealRelativeTorsion(const AtomKey &key);
	void initCrdPosi0(const std::map<std::string,NSPgeometry::XYZ> & crd0=std::map<std::string,NSPgeometry::XYZ>());
	void setBackBoneTorsionAt(unsigned int posi,const BackBoneTorsions & tor);
	void setBackBoneTorsions(const std::map<unsigned int,BackBoneTorsions> &tors);
	void setPhi(unsigned int posi, double phi,bool shift=false){
		assert(posi <length());
		AtomKey c_k=topo()->genKey(posi,BackBoneData::atomnames[BackBoneData::C]);
		AtomKey ca_k=topo()->genKey(posi,BackBoneData::atomnames[BackBoneData::CA]);
		if(!shift)
			phi =phi-intcrdmap_.at(c_k).torsion;
		rotateBond(ca_k,phi);
	}
	void setPsi(unsigned int posi, double psi,bool shift=false){
		assert(posi <length());
		AtomKey o_k=topo()->genKey(posi,BackBoneData::atomnames[BackBoneData::O]);
		AtomKey c_k=topo()->genKey(posi,BackBoneData::atomnames[BackBoneData::C]);
		psi = psi + 3.14159265358979323846;
		if(!shift) psi =psi-intcrdmap_.at(o_k).torsion;
		rotateBond(c_k,psi);
	}
	void setOmega(unsigned int posi, double omega,bool shift=false){
		assert(posi <length());
		AtomKey ca_k=topo()->genKey(posi,BackBoneData::atomnames[BackBoneData::CA]);
		AtomKey n_k=topo()->genKey(posi,BackBoneData::atomnames[BackBoneData::N]);
		if(!shift) omega =omega -intcrdmap().at(ca_k).torsion;
		rotateBond(n_k,omega);
	}

	NSPgeometry::XYZ getCrd(unsigned int posi, const std::string &atomname) const {
		assert(posi < length());
		AtomKey k=topo()->genKey(posi,atomname);
		return crdmap_.at(k);
	}
	BackBoneTorsions getBackBoneTorsions(unsigned int posi) const {
		return BackBoneTorsions(getPhi(posi),getPsi(posi),getOmega(posi));
	}
	double getPhi(unsigned int posi)const {
		assert(posi >= 0);
		AtomKey k=topo()->genKey(posi,BackBoneData::atomnames[BackBoneData::C]);
		return intcrdmap_.at(k).torsion;
	}
	double getPsi(unsigned int posi) const {
		assert(posi <length());
		if(posi <length()-1) {
			AtomKey k=topo()->genKey(posi+1,BackBoneData::atomnames[BackBoneData::N]);
			return intcrdmap_.at(k).torsion;
		} else {
			AtomKey k=topo()->genKey(posi,BackBoneData::atomnames[BackBoneData::O]);
			return intcrdmap_.at(k).torsion+3.14159265358979323846;
		}
	}
	double getOmega(unsigned int posi)const {
		assert(posi <length());
		AtomKey k=topo()->genKey(posi,BackBoneData::atomnames[BackBoneData::CA]);
		return intcrdmap().at(k).torsion;
	}
	NSPgeometry::XYZ getCrd(unsigned int posi, const std::string &atomname) {
		assert(posi <length());
		AtomKey k=topo()->genKey(posi,atomname);
		assert(crdValid(k));
		return crdmap_[k];
	}
	void getBackBoneCrd(unsigned int posi,std::vector<NSPgeometry::XYZ> &result){
		result.push_back(getCrd(posi,BackBoneData::atomnames[BackBoneData::N]));
		result.push_back(getCrd(posi,BackBoneData::atomnames[BackBoneData::CA]));
		result.push_back(getCrd(posi,BackBoneData::atomnames[BackBoneData::C]));
		result.push_back(getCrd(posi,BackBoneData::atomnames[BackBoneData::O]));
	}

/*	void getCrdCA2CA(unsigned int beginposi, unsigned int endposi,std::vector<NSPgeometry::XYZ> &result){
		for(unsigned int posi=beginposi;posi<=endposi;++posi) {
			if(posi != beginposi)result.push_back(getCrd(posi,BackBoneData::atomnames[BackBoneData::N]));
			result.push_back(getCrd(posi,BackBoneData::atomnames[BackBoneData::CA]));
			if(posi !=endposi) {
				result.push_back(getCrd(posi,BackBoneData::atomnames[BackBoneData::C]));
				result.push_back(getCrd(posi,BackBoneData::atomnames[BackBoneData::O]));
			}
		}
	}*/
	unsigned int length() const {return topo()->length();}
	void writePDB(std::ostream & os) const;
	void writeIntCrd(std::ostream &os) const ;
private:
};

}

#endif /* PROTEINREP_CHAINTREE_H_ */
