/*
 * backboneenergy.h
 *
 *  Created on: 2017年8月2日
 *      Author: hyliu
 */

#ifndef BACKBONE_BACKBONEENERGY_H_
#define BACKBONE_BACKBONEENERGY_H_
#include "backbone/backbonesite.h"
#include "backbone/backbonemoves.h"
#include "dataio/parameters.h"
#include "pdbstatistics/phipsidistr.h"
#include <memory>
namespace NSPproteinrep {
class BackBoneEnergy;
class ChainPack;
class ChainPackMoves;
typedef NSPdataio::TypedMultiInstanceControls<BackBoneEnergy> EnergyControls;


class BackBoneWindow {
public:
	BackBoneWindow(int ww) :
			windowwidth_(ww) {
		;
	}
	int windowwidth() const {
		return windowwidth_;
	}
	template<typename T>
	static int maxwindowwidth(const std::vector<T> &terms) {
		int res = -1;
		for (auto &t : terms) {
			if (res < t->windowwidth())
				res = t->windowwidth();
		}
		return res;
	}
protected:
	int windowwidth_;
};
struct EnergyComponents{
	enum CompTypes {TOTAL,PHIPSI,BLOCKLOCAL,CLASH,BLOCKPACKING,SPECIAL};
	static int NUM_ECOMP;
	std::vector<double> energies;
	EnergyComponents(){energies.resize(NUM_ECOMP,0.0);}
	void reset(){energies.clear();energies.resize(NUM_ECOMP,0.0);}
};

class ChainEnergyControl {
public:
	ChainEnergyControl(const std::string &parametersetname=std::string())
		:parasetname_(parametersetname),weights_(EnergyComponents::NUM_ECOMP,1.0){;}
	bool positionmasked(int posi) {if(positionmask_.empty()) return false; else return positionmask_[posi];}
	double & energy(int i) {return ecomp_.energies[i];}
	const double &energy(int i) const {return ecomp_.energies[i];}
	double eweight(int comptype) {return weights_[comptype];}
	std::vector<double> &eweights(){return weights_;}
	const std::vector<double>&eweights() const {return weights_;}
	std::string &origsscode() {return origsscode_;}
	const std::string &origsscode() const {return origsscode();}
	std::string & refpbseq() {return refpbseq_;}
	const std::string &refpbseq() const {return refpbseq_;}
	EnergyComponents &ecomp(){return ecomp_;}
	void setpositionmask(const std::vector<BackBoneSite> &chain);
	void resetecomp() {ecomp_.reset();}
	bool & ssaspbtype(){return ssaspbtype_;}
	bool &phipsi_ignoreresname(){return phipsi_ignoreresname_;}
	bool phipsi_ignoreresname() const {return phipsi_ignoreresname_;}
//	int packingminsep()const{ return blockpackingsep_;}
private:
	std::vector<bool> positionmask_;
	bool phipsi_ignoreresname_{false};
	EnergyComponents ecomp_;
	std::string refpbseq_;
	std::string origsscode_;
	std::vector<double> weights_;
	bool ssaspbtype_{false};
//	int blockpackingsep_{-1};
	std::string parasetname_;
};
int initenergycontrols(const std::string &ecid,const std::vector<std::string> &controlines=std::vector<std::string>());
inline int initenergycontrols(const std::vector<std::string> &controlines=std::vector<std::string>()){
	return initenergycontrols(std::string(),controlines);
}
int adjustenergycontrols(const std::string &ecid,const std::vector<std::string> &controlines);
inline int adjustenergycontrols(const std::vector<std::string> &controlines){
	return adjustenergycontrols(std::string(),controlines);
}
ChainEnergyControl prepareenergycontrol(std::vector<BackBoneSite> *chain,const std::string &parasetname);
class BackBoneLocalTerm: public BackBoneWindow {
public:

	BackBoneLocalTerm() :
			BackBoneWindow(1) {
	}
	BackBoneLocalTerm(int ww) :
			BackBoneWindow(ww) {
	}
	virtual double energy(
			std::vector<BackBoneSite>::const_iterator iter) const=0;
//	virtual double energy(std::vector<const BackBoneSite *> & window) const;
	virtual ~BackBoneLocalTerm() {
		;
	}
	std::string * &refpbseq() { return refpbseq_;}
	std::string * const &refpbseq() const {return refpbseq_;}
protected:
	std::string *refpbseq_{nullptr};
};

class PhiPsiTerm: public BackBoneLocalTerm {
public:
	PhiPsiTerm() :
			BackBoneLocalTerm(1) {
		;
	}
	PhiPsiTerm(bool ignore_resname): BackBoneLocalTerm(1),ignoreresname_(ignore_resname){;}
	virtual double energy(
			std::vector<BackBoneSite>::const_iterator iter) const {
//		std::cout <<"Phi-Psi for residue " <<iter->resseq <<std::endl;
		if(!refpbseq_->empty()){
			char expectedpbtype=(*refpbseq_)[iter->resseq];
			if(expectedpbtype =='m'){
				return NSPpdbstatistics::PhiPsiDistr::helixdistr().statisticalenergy(
						iter->phi(), iter->psi());
			} else if(expectedpbtype=='d'){
				return NSPpdbstatistics::PhiPsiDistr::stranddistr().statisticalenergy(
						iter->phi(), iter->psi());
			}
		}
		if(iter->sscode=='H') return NSPpdbstatistics::PhiPsiDistr::helixdistr().statisticalenergy(
				iter->phi(), iter->psi());
		if(iter->sscode=='E') return NSPpdbstatistics::PhiPsiDistr::stranddistr().statisticalenergy(
				iter->phi(), iter->psi());
		return energy(iter->phi(), iter->psi(), iter->resname);
	}
	double energy(double phi, double psi, const std::string &resname) const;
	virtual ~PhiPsiTerm() {
		;
	}
	bool ignoreresname_{false};
};

class PBlockTerm: public BackBoneLocalTerm {
public:
	PBlockTerm() :
			BackBoneLocalTerm(5) {;}
	virtual double energy(std::vector<BackBoneSite>::const_iterator iter) const;
	double energy(int resseq,char pbtype,const std::vector<double> &pbtorsions) const;
	virtual ~PBlockTerm() {
		;
	}
};

class BackBonePackingTerm: public BackBoneWindow {
public:
	BackBonePackingTerm(double cacut2, int minsep, int ww) :
			cacutoff2_(cacut2), minsep_(minsep),BackBoneWindow(ww) {;}
	virtual double energy(std::vector<BackBoneSite>::const_iterator iter1,
			std::vector<BackBoneSite>::const_iterator iter2, double *cadist2) const=0;
	virtual double energy(const BackBoneSite &s1,
			const BackBoneSite &s2, double *cadist2) const=0;
	virtual ~BackBonePackingTerm() {
		;
	}
protected:
	double cacutoff2_;
	int minsep_;
};
class StericClashTerm: public BackBonePackingTerm {
public:
	StericClashTerm() :
		BackBonePackingTerm(81.0, 2, 1) {
		auto & par=EnergyControls::getparameterset();
		if(!par.initialized) initenergycontrols();
		if(par.keydefined("ClashEnergy")) par.getval("ClashEnergy",&clashenergy_);
	}
	virtual double energy(std::vector<BackBoneSite>::const_iterator iter1,
			std::vector<BackBoneSite>::const_iterator iter2,double *cadist2) const;
	virtual double energy(const BackBoneSite &s1,
			const BackBoneSite &s2, double *cadist2) const;
	virtual ~StericClashTerm() {
		;
	}
	double clashenergy() const{return clashenergy_;}
private:
	double clashenergy_{10000.0};
};
class BlockPackingTerm: public BackBonePackingTerm {
public:
	BlockPackingTerm(int minsep) :
			BackBonePackingTerm(196.,minsep,5) {
		;
	}
	virtual double energy(std::vector<BackBoneSite>::const_iterator iter1,
			std::vector<BackBoneSite>::const_iterator iter2,double *cadist2) const;
	virtual double energy(const BackBoneSite &s1,
			const BackBoneSite &s2, double *cadist2) const;
	virtual ~BlockPackingTerm() {
		;
	}
};

class LocalBackBoneEnergy {
public:
	double totalenergy(const std::vector<BackBoneSite> & chain,ChainEnergyControl *ce,int ignorehead=0,
			int ignoretail=0);
	double partialE(const std::vector<BackBoneSite> &chain, int movestart,
			int moveend, const std::vector<BackBoneSite> &moved,ChainEnergyControl *ce);
	double partialE(const std::vector<BackBoneSite> &chain, int partstart,
			int partend,ChainEnergyControl *ce);
	double totalenergy(const ChainPack & cpk, EnergyComponents *ecomp);
	void partialE(const ChainPack &cpk, ChainPackMoves *moves,int mvidx=0);
	void addterm(std::shared_ptr<BackBoneLocalTerm> & newterm, int comptype) {
		terms_.push_back(newterm);
		comptypes_.push_back(comptype);
		wwmax_ = BackBoneWindow::maxwindowwidth(terms_);
		wwh_ = wwmax_ / 2;
	}
	void setrefpbseq(std::string *refpbseq){
		for (auto &t:terms_) {
			t->refpbseq()=refpbseq;
		}
	}
private:
	std::vector<std::shared_ptr<BackBoneLocalTerm> > terms_;
	std::vector<int> comptypes_;
	int wwmax_ { 0 };
	int wwh_ { 0 };
};

class BackBonePackingEnergy {
public:
	double totalenergy(const std::vector<BackBoneSite> & chain,ChainEnergyControl *ce,int ignorehead=0,
			int ignoretail=0);
	double totalenergy(const ChainPack & cpk, EnergyComponents *ecomp);
	void partialE(const ChainPack &cpk, ChainPackMoves *moves,int mvidx=0);
	double partialE(const std::vector<BackBoneSite> &chain, int movestart,
			int moveend, const std::vector<BackBoneSite> &moved,ChainEnergyControl *ce);
	double partialE(const std::vector<BackBoneSite> &chain, int partstart,
			int partend,ChainEnergyControl *ce);
	void addterm(std::shared_ptr<BackBonePackingTerm> & newterm, int comptype) {
		terms_.push_back(newterm);
		comptypes_.push_back(comptype);
		wwmax_ = BackBoneWindow::maxwindowwidth(terms_);
		wwh_ = wwmax_ / 2;
	}
private:
	std::vector<std::shared_ptr<BackBonePackingTerm> > terms_;
	std::vector<int> comptypes_;
	int wwmax_ { 0 };
	int wwh_ { 0 };
};
LocalBackBoneEnergy makelocalbackboneenergy(const std::vector<double> &weights
		=std::vector<double>(EnergyComponents::NUM_ECOMP,1.0),bool phipsi_ignoreresname=false);
BackBonePackingEnergy makebackbonepackingenergy(
		const std::vector<double> &weights
				=std::vector<double>(EnergyComponents::NUM_ECOMP,1.0));
class BackBoneEnergy{
public:

	BackBoneEnergy(){
		lene_=makelocalbackboneenergy();
		pene_= makebackbonepackingenergy();
	}
	BackBoneEnergy(const ChainEnergyControl &ce){
		lene_=makelocalbackboneenergy(ce.eweights(),ce.phipsi_ignoreresname());
		pene_=makebackbonepackingenergy(ce.eweights());
	}
	double totalenergy(const std::vector<BackBoneSite> & chain,ChainEnergyControl *ce){
		ce->resetecomp();
		lene_.setrefpbseq(&ce->refpbseq());
		double etot=lene_.totalenergy(chain,ce)+pene_.totalenergy(chain,ce);
		ce->energy(EnergyComponents::TOTAL)=etot;
		return etot;
	}
	double totalenergy(const ChainPack & cpck,EnergyComponents *ecomp);
	double deltaE(const ChainPack & cpck,ChainPackMoves *moves, int mvidx=0);
	double totalenergy(const std::vector<BackBoneSite> & chain){
		ChainEnergyControl ce;
		return totalenergy(chain,&ce);
	}
	double deltaE(BackBoneMoves &moves,int mvidx,ChainEnergyControl *ce,
			EnergyComponents *ecompa, EnergyComponents *ecompb);
	double deltaE(BackBoneMoves &moves,int mvidx){
		ChainEnergyControl ce;
		EnergyComponents ecompa,ecompb;
		return deltaE(moves,mvidx,&ce,&ecompa,&ecompb);
	}
	double partialE(const std::vector<BackBoneSite> &chain, int partstart,
			int partend,ChainEnergyControl *ce);
private:
	LocalBackBoneEnergy lene_;
	BackBonePackingEnergy pene_;
};
}

#endif /* BACKBONE_BACKBONEENERGY_H_ */
