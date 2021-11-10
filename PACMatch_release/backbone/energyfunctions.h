/*
 * energyfunctions.h
 *
 *  Created on: 2017年4月25日
 *      Author: hyliu
 */

#ifndef BACKBONE_ENERGYFUNCTIONS_H_
#define BACKBONE_ENERGYFUNCTIONS_H_

#include "backbone/backbonesite.h"
#include <map>
namespace NSPproteinrep {

template<typename ENETERMTYPE>
ENETERMTYPE & eneterminstance(double w=1.0) {
	static ENETERMTYPE ene;
	ene.weight()=w;
	return ene;
}
class WeightedTerm {
public:
	virtual std::string termname()=0;
	double &weight() {return weight_;}
	const double &weight() const {return weight_;}
	virtual ~WeightedTerm(){;}
private:
	double weight_{1.0};
};
class SiteEne:public WeightedTerm {
public:
	SiteEne() {
		;
	}
	virtual double energy(const BackBoneSite & bs)=0;
	virtual ~SiteEne() {
		;
	}
};

class SitePairEne:public WeightedTerm {
public:
	SitePairEne() {
		;
	}
	virtual double energy(const BackBoneSite & bs1, const BackBoneSite &bs2,
			int sep = 10000)=0;
	virtual ~SitePairEne() {
		;
	}
};

class StericClash: public SitePairEne {
public:
	friend StericClash & eneterminstance<StericClash>(double w);
	virtual double energy(const BackBoneSite & bs1, const BackBoneSite &bs2,
			int sep = 10000) {
		if (sep == 1){
			int clash=nextresidue_atomsclashed(bs1,bs2);
			if(clash >=2) return eclash_;
			else if(clash==1) return enrclash_;
			return 0.0;
		}
		if (atomsclashed(bs1, bs2))
			return eclash_;
		return 0.0;
	}
	double eclash() const {return eclash_;}
	void setenrclash(double e) {enrclash_=e;}
	virtual std::string termname(){return "StericClash";}
protected:
	StericClash() :
			minsep_(2), eclash_ ( 1.e20 ),enrclash_(100.0) {
		;
	}
private:
	int minsep_;
	double eclash_;
	double enrclash_;
};

class SiteAtomContact: public SitePairEne {
public:

	static SiteAtomContact & instance(const std::string &atom1, const std::string & atom2,
			double rcut,double w=1.0) {
		static SiteAtomContact inst;
		inst.init(atom1, atom2,rcut);
		inst.weight()=w;
		return inst;
	}
	virtual std::string termname(){return "SiteAtomContact";}
	virtual double energy(const BackBoneSite & bs1, const BackBoneSite &bs2,
			int sep = 10000) ;
private:
	SiteAtomContact():minsep_(2) {
		;
	}
	void init(const std::string & atom1,const std::string &atom2,double rcut) {
		std::map<std::string,int> atomids{{"N",0},{"CA",1},{"C",2},{"O",3}};
		atom1_=atomids.at(atom1);
		atom2_=atomids.at(atom2);
		rcut2_ = rcut * rcut;
	}

private:
	int minsep_ {2};
	int atom1_{-1};
	int atom2_{-1};
	double rcut2_ { 0.0 };
};

class LocalEne:public WeightedTerm {
public:
	static int fragmentlength;
	LocalEne() {
		;
	}
	virtual double redundancy(int posi)=0;
	virtual double totalredundancy()=0;
	virtual int leftextension()=0;
	virtual int rightextension()=0;
	virtual bool defined(const std::vector<const BackBoneSite*> & sites)=0;
	virtual double energy(const std::vector<const BackBoneSite*> & sites)=0;
	virtual ~LocalEne() {
		;
	}
};
class TorsionEne: public LocalEne {
public:
	virtual double energy(const std::vector<const BackBoneSite *> & sites);
	virtual ~TorsionEne() {
		;
	}
	virtual bool defined(const std::vector<const BackBoneSite*> & sites){
		if(sites[fragmentlength/2]) return true;
		return false;
	}
	virtual std::string termname(){return "TorsionEne";}
	virtual double redundancy(int posi) {return redundancy_[posi];}
	virtual double totalredundancy() {return 1.0;}
	virtual int leftextension() {return 0;}
	virtual int rightextension() {return 0;}
	friend TorsionEne & eneterminstance<TorsionEne>(double w);
protected:
	TorsionEne();
private:
	std::vector<double> redundancy_{0.0,0.0,1.0,0.0,0.0};
};

class TorsionVecEne: public LocalEne {
public:
	virtual double energy(const std::vector<const BackBoneSite *> & sites);
	virtual ~TorsionVecEne() {
		;
	}
	virtual bool defined(const std::vector<const BackBoneSite*> & sites){
		for(int i=1; i<=motiflength_;++i) {
			if (sites[i]==nullptr) return false;
		}
		return true;
	}
	virtual std::string termname(){return "TorsionVecEne";}
	virtual double redundancy(int posi) {return redundancy_[posi];}
	virtual double totalredundancy() {return 4.0;}
	virtual int leftextension() {return 2;}
	virtual int rightextension() {return 1;}
	friend TorsionVecEne & eneterminstance<TorsionVecEne>(double w);
protected:
	TorsionVecEne();
private:
	std::vector<double> redundancy_{0.0,1.0,1.0,1.0,1.0};
	int motiflength_{4};
};

class PackingEne:public WeightedTerm {
public:
	PackingEne() {
		;
	}
	virtual double energy(const std::vector<const BackBoneSite *> &frag1,
			const std::vector<const BackBoneSite *> &frag2, int seqsep)=0;
	virtual ~PackingEne() {
		;
	}
};

class TetraSefEne: public PackingEne {
public:
	virtual double energy(const std::vector<const BackBoneSite *> &frag1,
			const std::vector<const BackBoneSite *> &frag2, int seqsep = 10000);
	virtual ~TetraSefEne() {
		;
	}
	virtual std::string termname(){return "TetraSefEne";}
	friend TetraSefEne & eneterminstance<TetraSefEne>(double w);
protected:
	TetraSefEne();
	int minsep_;
};

class EnergyTerms {
public:
	std::vector<SiteEne *> siteene_terms;
	std::vector<SitePairEne *> sitepairene_terms;
	std::vector<LocalEne *> localene_terms;
	std::vector<PackingEne *> packingene_terms;
	void addstericclash(double enrclash=100.0){
		sitepairene_terms.push_back(&eneterminstance<StericClash>());
		auto p=(StericClash *) (sitepairene_terms.back());
		p->setenrclash(enrclash);
	}
	void addtorsionene(double w=1.0) {
		localene_terms.push_back(&eneterminstance<TorsionEne>(w));
	}
	void addtorsionvecene(double w=1.0) {
		localene_terms.push_back(&eneterminstance<TorsionVecEne>(w));
	}
	void addtetrasefene(double w=1.0) {
		packingene_terms.push_back(&eneterminstance<TetraSefEne>(w));
	}
	void setemaxallowed(double e) {emaxallowed=e;}
	double emaxallowed{1.e15};
};
class MainChain;
double calctotalenergy(MainChain &mc, EnergyTerms &terms,
		std::map<std::string, double> & energies,double egiveup);
double calctotalenergy(MainChain &mc,
		const std::vector<SiteEne *> & terms,
		std::map<std::string, double> & energies,double egiveup);
double calctotalenergy(MainChain &mc,
		const std::vector<SitePairEne *> & terms,
		std::map<std::string, double> & energies,double egiveup);
double calctotalenergy(MainChain &mc,
		const std::vector<LocalEne *> & localterms,
		const std::vector<PackingEne *> & packingterms,
		std::map<std::string, double> & energies,double egiveup);
double calcloopenergy(MainChain &mc, int loopstart, int looplength,
		const std::vector<BackBoneSite> &newloop, EnergyTerms &terms,
		std::map<std::string, double> & energies,double egiveup);
double calcloopenergy(MainChain &mc,int loopstart, int looplength,
		const std::vector<BackBoneSite> &newloop,
		const std::vector<SiteEne *> & terms,
		std::map<std::string, double> & energies,double egiveup);
double calcloopenergy(MainChain &mc,int loopstart, int looplength,
		const std::vector<BackBoneSite> &newloop,
		const std::vector<SitePairEne *> & terms,
		std::map<std::string, double> & energies,double egiveup);
double calcloopenergy(MainChain &mc,int loopstart, int looplength,
		const std::vector<BackBoneSite> &newloop,
		const std::vector<LocalEne *> & localterms,
		const std::vector<PackingEne *> & packingterms,
		std::map<std::string, double> & energies,double egiveup);
double looploopenergy(const std::vector<BackBoneSite> &loop1,
		const std::vector<BackBoneSite> &loop2,
		SitePairEne * eterm,double egiveup);
void decomposeenergy(MainChain &mc,const EnergyTerms &terms,std::ostream &os);
}

#endif /* BACKBONE_ENERGYFUNCTIONS_H_ */
