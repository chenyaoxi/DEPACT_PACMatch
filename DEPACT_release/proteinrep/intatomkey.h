/*
 * intatomkey.h
 *
 *  Created on: 2016年11月15日
 *      Author: hyliu
 */

#ifndef PROTEINREP_INTATOMKEY_H_
#define PROTEINREP_INTATOMKEY_H_
#include <string>
#include <map>
#include <set>
#include <climits>
#include <cassert>
namespace NSPproteinrep {

extern const std::map<std::string,unsigned int> atomcodes;
extern const std::set<std::string> mainchainatoms;
extern const std::set<std::string> hbdonoratoms;
extern const std::set<std::string> hbacceptoratoms;
extern const std::map<char,std::string> residuecharnames;
extern const std::set<std::string> positiveresidues;
extern const std::set<std::string> negativeresidues;
extern const std::set<std::string> polarresidues;
extern const std::set<std::string> aromaticresidues;
extern const std::map<std::string, unsigned int> residuenamecodes;
typedef unsigned int AtomKeyType;

template <unsigned int ROTDIM, unsigned int RESDIM, unsigned int ATOMNMDIM,
		  unsigned int POSIDIM>

class IntAtomKey {
public:
	typedef AtomKeyType Key;
	static Key genKey(unsigned int posinumber, const std::string & atomname, unsigned int chainnum=0,
			const std::string & restype="ANY", unsigned int rotnum=0){
		auto it=atomcodes.find(atomname);
		assert (it != atomcodes.end());
		auto itr=residuenamecodes.find(restype);
		assert(itr != residuenamecodes.end());
		assert(posinumber < POSIDIM);
		assert(rotnum < ROTDIM);
		assert(ROTDIM*RESDIM*ATOMNMDIM*POSIDIM < UINT_MAX/(chainnum+1U));
		return rotnum+ROTDIM*(itr->second + RESDIM*(it->second+ATOMNMDIM*(posinumber + POSIDIM*chainnum)));
	}
	static unsigned int chainNumber(Key key){
		return key/(ROTDIM*RESDIM*ATOMNMDIM*POSIDIM);
	}
	static unsigned int posiNumber(Key key){
			return (key%(ROTDIM*RESDIM*ATOMNMDIM*POSIDIM))/(ROTDIM*RESDIM*ATOMNMDIM);
	}
	static unsigned int atomNmNumber(Key key){
			return (key%(ROTDIM*RESDIM*ATOMNMDIM))/(ROTDIM*RESDIM);
	}
	static unsigned int resnameCode(Key key){
			return (key%(ROTDIM*RESDIM))/ROTDIM;
	}
	static unsigned int rotCodeNumber(Key key){
			return key%ROTDIM;
	}
	static std::string atomName(Key k) {
		k= atomNmNumber(k);
		for(auto &nm:atomcodes)	if( nm.second == k) return nm.first;
	}
	static std::string residueName(Key k){
		k=resnameCode(k);
		for(auto &c:residuenamecodes) if(k==c.second) return c.first;
	}
	static bool positivelyCharged(const std::string & nm){
		return (positiveresidues.find(nm) != positiveresidues.end());
	}
	static bool positivelyCharged(Key k) {
		std::string nm=residueName(k);
		return positivelyCharged(nm);
	}
	static bool negativelyCharged(const std::string & nm){
		return (negativeresidues.find(nm) != negativeresidues.end());
	}
	static bool negativelyCharged(Key k) {
		std::string nm=residueName(k);
		return negativelyCharged(nm);
	}
	static bool isCharged(const std::string &nm) {
			return negativelyCharged(nm) || positivelyCharged(nm);
	}
	static bool isCharged(Key k) {
		std::string nm=residueName(k);
		return isCharged(nm);
	}
	static bool isPolar(const std::string &nm) {
		return isCharged(nm) || (polarresidues.find(nm) != polarresidues.end());
	}
	static bool isPolar(Key k) {
		std::string nm=residueName(k);
		return isPolar(nm);
	}
	static bool isAromatic(const std::string &nm) {
		return aromaticresidues.find(nm) != aromaticresidues.end();
	}
	static bool isAromatic(Key k) {
		std::string nm=residueName(k);
		return isAromatic(nm);
	}
	static bool isMainChain(const std::string &nm) {
		return mainchainatoms.find(nm) != mainchainatoms.end();
	}
	static bool isMainChain(Key k) {
		std::string nm=atomName(k);
		return isMainChain(nm);
	}
	static bool isHBDonor(const std::string &nm) {
		return hbdonoratoms.find(nm) != hbdonoratoms.end();
	}
	static bool isHBDonor(Key k) {
		std::string nm=atomName(k);
		return isHBDonor(nm);
	}
	static bool isHBAcceptor(const std::string &nm) {
		return hbacceptoratoms.find(nm) != hbacceptoratoms.end();
	}
	static bool isHBAcceptor(Key k) {
		std::string nm=atomName(k);
		return isHBAcceptor(nm);
	}
};
typedef IntAtomKey<220,50,80,1000> AtomKeyTypeR; //rotamer based key
typedef IntAtomKey<1,50,80,4000> AtomKeyTypeA; //residue based key
}


#endif /* PROTEINREP_INTATOMKEY_H_ */
