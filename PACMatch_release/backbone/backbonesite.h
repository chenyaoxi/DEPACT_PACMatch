/*
 * backbonesite.h
 *
 *  Created on: 2016年3月9日
 *      Author: hyliu
 */

#ifndef BACKBONESITE_H_
#define BACKBONESITE_H_
#include "geometry/calculators.h"
#include "geometry/rotation.h"
#include "dataio/inputlines.h"
//#include "pdbstructure.h"
//#include "chainattrib.h"
#include "geometry/quaternioncrd.h"
#include "geometry/localframe.h"
#include "geometry/relativeposition.h"
#include "proteinrep/pdbrecord.h"
//#include "hbgeometry.h"
#include <iostream>
#include <fstream>
#include <map>

#include <boost/algorithm/string.hpp>
namespace NSPproteinrep {
struct BackBoneSite {
	enum DataField {
		PHI = 0,
		PSI = 1,
		OMIGA = 2,
		SASA = 3,
		NCRD = 4,
		CACRD = 7,
		CCRD = 10,
		OCRD = 13,
		DATADIM = 17
	};
	typedef std::vector<BackBoneSite> backbonesite_attrib_type;
//	static std::string attrib_name;

	double *data() {
		return data_;
	}
	double phi() const {
		return data_[PHI];
	}
	double phi(const BackBoneSite &pevioussite);
	double psi(const BackBoneSite &nextsite);
	double omiga(const BackBoneSite &nextsite);
	double psi() const {
		return data_[PSI];
	}
	double omiga() const {
		return data_[OMIGA];
	}
	bool nextpepcis() const {
		return(data_[OMIGA]>-90.0 && data_[OMIGA]<90.0);
	}
	void settorsions(double phi,double psi){
		data_[PHI]=phi;data_[PSI]=psi;
		newocrdfrompsi();
	}
	const double *data() const {
		return data_;
	}
	char sscodechar() const;
	std::string toString() const;
	NSPgeometry::LocalFrame localframe();
	template<typename IT>
	static BackBoneSite read(IT &lineit) {
		BackBoneSite bs;
		std::vector<std::string> line1 = NSPdataio::parseline(*lineit++,
				std::vector<int>());
		bs.pdbid = boost::trim_copy(line1[0]);
		bs.chainid = line1[1][0];
		bs.resid = boost::lexical_cast<int>(line1[2]);
		bs.resname = boost::trim_copy(line1[3]);
		bs.resseq = boost::lexical_cast<int>(line1[4]);
		bs.sscode = line1[5][0];
		for (int i = 0; i < NCRD; ++i) {
			bs.data_[i] = boost::lexical_cast<double>(line1[6 + i]);
		}
		for (int i = NCRD; i < NCRD + 12; i += 3) {
			line1 = NSPdataio::parseline(*lineit++, std::vector<int>());
			;
			for (int j = 0; j < 3; j++)
				bs.data_[i + j] = boost::lexical_cast<double>(line1[j]);
		}
		return bs;
	}
	NSPgeometry::XYZ getcrd(int i) const {
		return NSPgeometry::XYZ(data_[i], data_[i + 1], data_[i + 2]);
	}
	NSPgeometry::XYZ ncrd() const {
		return getcrd(NCRD);
	}
	NSPgeometry::XYZ cacrd() const {
		return getcrd(CACRD);
	}
	NSPgeometry::XYZ ccrd() const {
		return getcrd(CCRD);
	}
	NSPgeometry::XYZ ocrd() const {
		return getcrd(OCRD);
	}
	NSPgeometry::XYZ changecrd(const std::vector<double> &crd) {
		assert(crd.size() >= 12);
		for (int i = 0; i < 12; ++i)
			data_[NCRD + i] = crd[i];
	}
	NSPgeometry::XYZ hcrd() const {
		NSPgeometry::XYZ c=getcrd(CCRD);
		NSPgeometry::XYZ ca=getcrd(CACRD);
		NSPgeometry::XYZ n=getcrd(NCRD);
		double theta=120.0*3.14359265/180.0;
		double t=phi()*3.14159265/180.0+3.14159265;
		return NSPgeometry::InternaltoXYZ(n,ca,c,1.0,theta,t);
	}
	NSPgeometry::XYZ cd_procrd() const {
		NSPgeometry::XYZ c=getcrd(CCRD);
		NSPgeometry::XYZ ca=getcrd(CACRD);
		NSPgeometry::XYZ n=getcrd(NCRD);
		double theta=120.0*3.14359265/180.0;
		double t=phi()*3.14159265/180.0+3.14159265;
		return NSPgeometry::InternaltoXYZ(n,ca,c,1.47,theta,t);
	}
	NSPgeometry::XYZ cbcrd(double b0=1.5) const {
		NSPgeometry::XYZ c=getcrd(CCRD);
		NSPgeometry::XYZ ca=getcrd(CACRD);
		NSPgeometry::XYZ n=getcrd(NCRD);
		double theta=109.5*3.14159265/180.0;
		double t=120*3.14159265/180.0;
		return NSPgeometry::InternaltoXYZ(ca,c,n,b0,theta,t);
	}
	NSPgeometry::XYZ changecrd(const std::vector<NSPgeometry::XYZ> &crd) {
		assert(crd.size() >= 4);
		int indx = NCRD;
		for (int i = 0; i < 4; i++) {
			data_[indx++] = crd[i].x_;
			data_[indx++] = crd[i].y_;
			data_[indx++] = crd[i].z_;
		}
	}
	void newocrdfrompsi(){
		NSPgeometry::Rotation r(ncrd(),cacrd(),ccrd(),ocrd(),(data_[PSI]+180.0)*3.14159265/180.0);
		NSPgeometry::XYZ newo=r.applytoCopy(ocrd());
		data_[OCRD]=newo.x_;
		data_[OCRD+1]=newo.y_;
		data_[OCRD+2]=newo.z_;
	}
	void newncrdfrompsi(){
		NSPgeometry::Rotation r(ocrd(),ccrd(),cacrd(),ncrd(),(data_[PSI]+180.0)*3.14159265/180.0);
		NSPgeometry::XYZ newn=r.applytoCopy(ncrd());
		data_[OCRD]=newn.x_;
		data_[OCRD+1]=newn.y_;
		data_[OCRD+2]=newn.z_;
	}
	template<typename ITER>
	void changecrd(ITER iter) {
		int indx = NCRD;
		for (int i = 0; i < 4; i++) {
			data_[indx++] = iter->x_;
			data_[indx++] = iter->y_;
			data_[indx++] = iter->z_;
			++iter;
		}
	}
	void translate(NSPgeometry::XYZ t){
		int index=NCRD;
		for (int i = 0; i < 4; ++i){
			data_[index++] += t.x_;
			data_[index++] += t.y_;
			data_[index++] += t.z_;
		}

	}
	void rotate(const NSPgeometry::Rotation & r) {
		std::vector<NSPgeometry::XYZ> crd;
		getcrd(crd);
		for(auto &p:crd) {
			r.apply(&p);
		}
		changecrd(crd.begin());
	}
	void getcrd(std::vector<double> & crd) const {
		for (int i = NCRD; i < NCRD + 12; ++i)
			crd.push_back(data_[i]);
	}
	void getcrd(std::vector<NSPgeometry::XYZ> & crd) const {
		for (int i = NCRD; i < NCRD + 12; i += 3)
			crd.push_back(NSPgeometry::XYZ(data_[i], data_[i + 1], data_[i + 2]));
	}

	void genPdbRecords(std::vector<PdbRecord> *records,int posi,int natmid=1) const;

	std::string pdbid;
	std::string resname;
	char chainid;
	char sscode;
	int resid;
	int resseq;
	bool isgap{false};
	double data_[DATADIM];
};

class BBAtomOrder{
public:
	enum {H=0,N=1,CA=2,C=3,O=4,CB=5,NATOMS=6};
};

class BBatomdis2 {

public:
	BBatomdis2(BackBoneSite &s1,BackBoneSite &s2);
	double operator()(int a1,int a2){return dist2_[a1*5+a2];}
	std::vector<double> & dist2() {return dist2_;}
private:
	std::vector<double> dist2_;
};


/*
void make_backbonesites_chain(ChainData & chain,
		ChainAttrib<BackBoneSite::backbonesite_attrib_type> &sites);
*/
template<typename SiteIterator>
bool sameSSsegment(SiteIterator it1, SiteIterator it2) {
	char sscode = it1->sscode;
	int resid_old = it1->resid;
	if (sscode != 'H' && it1->sscode != 'E')
		return false;
	do {
		it1++;
		if (it1->sscode != sscode || it1->resid - resid_old > 3)
			return false;
//		std::cout <<it1->sscode << resid_old <<" " <<it1->resid;
		resid_old = it1->resid;
	} while (it1 != it2);
	return true;
}
template<typename SiteIterator>
bool fragstartsite(SiteIterator iter, SiteIterator end, int length,
		const std::string & SSstring = std::string(),bool checkterminus=true) {
	const BackBoneSite *bs1 = &(*iter);
	if(checkterminus) {
		if (bs1->phi() == 360.0 || bs1->psi() == 360.0 || bs1->omiga() == 360.0)
			return false;
	}
	std::string pdbid = bs1->pdbid;
	if(checkterminus) {
		if (chainstartsite(iter) || chainendsite(iter, end))
			return false;  //ignore fragment containing first position
	}
	if (!SSstring.empty()) {
		assert(SSstring.length() == length);
		if (SSstring[0] != 'X') {
			if (SSstring[0] == 'H' || SSstring[0] == 'E') {
				if (bs1->sscode != SSstring[0])
					return false;
			} else if (bs1->sscode == 'H' || bs1->sscode == 'E')
				return false;
		}
	}
	char chainid = bs1->chainid;
	int resid = bs1->resid;
	for (int l = 0; l < length - 1; ++l) {
		bool broken = false;
		const BackBoneSite *bs0 = bs1;
		bs1 = &(*(++iter));
		if(checkterminus) {
			if (bs1->phi() == 360.0 || bs1->psi() == 360.0 || bs1->omiga() == 360.0)
				return false;
			if (chainendsite(iter, end))
				return false; //ignore fragment containing last position
		}
		if (!SSstring.empty()) {
			assert(SSstring.length() == length);
			if (SSstring[l + 1] != 'X') {
				if (SSstring[l + 1] == 'H' || SSstring[l + 1] == 'E') {
					if (bs1->sscode != SSstring[l + 1])
						return false;
				} else if (bs1->sscode == 'H' || bs1->sscode == 'E')
					return false;
			}
		}
		/*		if( iter == end || bs1->pdbid != pdbid
		 || bs1->chainid != chainid
		 ||bs1->resseq != resseq+l+1) {
		 //			std::cout <<pdbid <<bs1->pdbid <<chainid <<bs1->chainid
		 //					<<resseq<<bs1->resseq <<std::endl;
		 return false;
		 }
		 */
		const double *d0 = bs0->data();
		NSPgeometry::XYZ c0(d0[BackBoneSite::CCRD], d0[BackBoneSite::CCRD + 1],
				d0[BackBoneSite::CCRD + 2]);
		const double *d1 = bs1->data();
		NSPgeometry::XYZ n1(d1[BackBoneSite::NCRD], d1[BackBoneSite::NCRD + 1],
				d1[BackBoneSite::NCRD + 2]);
		double bond = NSPgeometry::distance(c0, n1);
		if (bond > 1.60) {
			//	if(bs1->resid != resid+l+1) {
//		std::cout <<l<<" "<<pdbid <<bs1->pdbid <<chainid <<bs1->chainid
//					<<resid<<" "<<bs1->resid <<std::endl;
//			std::cout <<c0.x_ <<"\t" <<c0.y_ <<"\t" <<c0.z_ <<"\n";
//			std::cout <<n1.x_ <<"\t" <<n1.y_ <<"\t" <<n1.z_ <<"\n";
//			std::cout <<bond <<std::endl;
			return false;
		}
	}
	return true;
}
template<typename SiteIterator>
bool chainendsite(SiteIterator iter, SiteIterator end) {
	if(iter == end | iter+1 == end) return true;
	const BackBoneSite* bs = &(*iter);
	const BackBoneSite* bs1 = &(*(++iter));
	if (bs1->pdbid != bs->pdbid || bs1->chainid != bs->chainid)
		return true;
	return false;
}
template<typename SiteIterator>
bool chainstartsite(SiteIterator iter) {
	const BackBoneSite *bs = &(*iter);
	const BackBoneSite *bs1 = &(*(iter - 1));

	if (bs1->pdbid != bs->pdbid || bs1->chainid != bs->chainid)
		return true;
	return false;
}
template<typename SiteIterator>
SiteIterator findchainstart(SiteIterator iter, SiteIterator begin) {
	while (iter != begin) {
		if (!chainstartsite(iter))
			--iter;
		else
			return iter;
	}
	return begin;
}

template<typename SiteIterator>
SiteIterator findnextchainstart(SiteIterator iter, SiteIterator end) {
	while (iter != end) {
		if (!chainendsite(iter, end))
			++iter;
		else
			return ++iter;
	}
	return end;
}
template<typename SiteIterator>
std::string getSSstring(SiteIterator begin,int length) {
	std::string SS(length,' ');
	for (int i = 0; i < length; ++i) {
		auto iter = begin + i;
		SS[i] = iter->sscode;
		if (SS[i] != 'H' && SS[i] != 'E')
			SS[i] = 'C';
	}
	return SS;
}
template<typename SiterIterator>
void gettorsionvector(SiterIterator begin, int length, std::vector<double> &d) {
	for (int i = 0; i < length; ++i) {
		auto iter = begin + i;
		double ang = iter->phi();
		if (ang > 180.0)
			ang -= 360.0;
		if (ang < -180.0)
			ang += 360.0;
		d.push_back(ang);
		ang = iter->psi();
		if (ang > 180.0)
			ang -= 360.0;
		if (ang < -180.0)
			ang += 360.0;
		d.push_back(ang);
		if (i < length - 1) {
			ang = iter->omiga();
			if (ang > 180.0)
				ang -= 360.0;
			if (ang < -180.0)
				ang += 360.0;
			d.push_back(ang);
		}
	}
}
void readbackbonesites(const std::string & filename,
		std::vector<BackBoneSite> &sites);
void writeSitesToPDB(std::ostream & os,const std::vector<BackBoneSite> & sites);
bool readbackbonesites(std::istream &is, int number, std::vector<BackBoneSite> &sites);
NSPgeometry::QuaternionCrd relativeorientation(BackBoneSite &s1, BackBoneSite &s2);
NSPgeometry::RelativePosition relativeposition(BackBoneSite &s1, BackBoneSite &s2);
void readbackbonefrompdb(const std::string &filename,
		std::vector<std::vector<BackBoneSite>> &chains);

std::vector<BackBoneSite> generaterandombackbone(int length,
		const std::vector<std::pair<int,int>> &helixregions,
		const std::vector<std::pair<int,int>> &strandregions);

std::vector<int> findcissites(const std::vector<BackBoneSite> &chain);
//BackBoneHB hbproperty(BackBoneSite &s1, BackBoneSite &s2);
bool atomsclashed(const BackBoneSite &s1, const BackBoneSite &s2);
int nextresidue_atomsclashed(const BackBoneSite &s1, const BackBoneSite &s2);
void genbackbonesite(BackBoneSite *psite, bool cispep,
		double phi, double psi, BackBoneSite *newsite);
void genprevbackbonesite(BackBoneSite *nsite,double omiga,double psi, double phi,
		BackBoneSite *newsite);
std::vector<BackBoneSite> loops2gaps(const std::vector<BackBoneSite> &chain,int extension=2);
std::vector<BackBoneSite> insertgaps(const std::vector<BackBoneSite> &chain,
		const std::map<int,int> &gapposi_lens);
std::vector<std::pair<int,int>> alignaspointsets(const std::vector<BackBoneSite> &seta,
		const std::vector<BackBoneSite> &setb,double rcut=3.0);
int sselements(const std::vector<BackBoneSite> &chain, std::vector<int> *elmtidseq,
		std::vector<std::pair<int,int>> *elmtposlens);
std::vector<double> extractcrd(const std::vector<BackBoneSite> &chain);
void assigncrd(const std::vector<double> &crd, std::vector<BackBoneSite> &chain);
}
#endif /* BACKBONESITE_H_ */
