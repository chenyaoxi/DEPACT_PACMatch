/*
 * backbonesite.cpp
 *
 *  Created on: 2016年3月9日
 *      Author: hyliu
 */
#include <backbone/backbonesite.h>
#include "geometry/calculators.h"
#include "geometry/localframe.h"

#include "geometry/calculators.h"
#include "proteinrep/idealgeometries.h"
#include "geometry/structalign.h"
#include "pdbstatistics/phipsidistr.h"
#include "dstl/randomengine.h"
#include "proteinrep/pdbreader.h"
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace NSPproteinrep;
//std::string BackBoneSite::attrib_name{"backbonesite"};
double BackBoneSite::phi(const BackBoneSite &ps) {
	double t = NSPgeometry::torsion(ps.getcrd(CCRD), getcrd(NCRD),
			getcrd(CACRD), getcrd(CCRD));
	double rad = 180.0 / 3.14159265358979323846;
	t *= rad;
	while (t > 180.0)
		t -= 360.0;
	while (t < -180.0)
		t += 360.0;
	data_[PHI] = t;
	return t;
}
double BackBoneSite::psi(const BackBoneSite &ns) {
	double t = NSPgeometry::torsion(getcrd(NCRD), getcrd(CACRD), getcrd(CCRD),
			ns.getcrd(NCRD));
	double rad = 180.0 / 3.14159265358979323846;
	t *= rad;
	while (t > 180.0)
		t -= 360.0;
	while (t < -180.0)
		t += 360.0;
	data_[PSI] = t;
	return t;
}
double BackBoneSite::omiga(const BackBoneSite &ns) {
	double t = NSPgeometry::torsion(getcrd(CACRD), getcrd(CCRD),
			ns.getcrd(NCRD), ns.getcrd(CACRD));
	double rad = 180.0 / 3.14159265358979323846;
	t *= rad;
	while (t > 180.0)
		t -= 360.0;
	while (t < -180.0)
		t += 360.0;
	data_[OMIGA] = t;
	return t;
}

std::string BackBoneSite::toString() const {
	std::ostringstream oss;
	oss << std::setw(5) << pdbid << std::setw(2) << chainid << std::setw(5)
			<< resid << std::setw(5) << resname << std::setw(5) << resseq
			<< std::setw(2) << sscode << " ";
	for (int i = 0; i < NCRD; ++i) {
		oss << std::setw(10) << std::setiosflags(std::ios::right)
				<< std::setiosflags(std::ios::fixed) << std::setprecision(2)
				<< data_[i];
	}
	std::string space(20, ' ');
	oss << std::endl;
	for (int i = NCRD; i < NCRD + 12; i += 3) {
		oss << space;
		for (int j = 0; j < 3; ++j)
			oss << std::setw(10) << std::setiosflags(std::ios::right)
					<< std::setiosflags(std::ios::fixed) << std::setprecision(3)
					<< data_[i + j];
		oss << std::endl;
	}
	return oss.str();
}
void BackBoneSite::genPdbRecords(std::vector<PdbRecord> *records, int posi,
		int natmid) const {
	std::vector<std::string> atomnames { "N", "CA", "C", "O" };
	for (int i = 0; i < 4; ++i) {
		records->push_back(PdbRecord());
		PdbRecord & record = records->back();
		record.label = "ATOM";
		record.chainid = chainid;
		record.atomname = atomnames[i];
		record.namesymbol = record.atomname.substr(0, 1);
		record.elementname[1] = record.namesymbol[0];
		record.namemodifier = record.atomname.substr(1);
		record.residuename = resname;
		record.atomid = natmid + i;
		record.residueid = posi;
		NSPgeometry::XYZ crd = getcrd(NCRD + 3 * i);
		record.x = crd.x_;
		record.y = crd.y_;
		record.z = crd.z_;
	}
}
void NSPproteinrep::writeSitesToPDB(std::ostream &os,
		const std::vector<BackBoneSite> & sites) {
	unsigned int nid;
	std::vector<PdbRecord> records;
	int posi = 0;
	for (auto & s : sites) {
		nid = records.size() + 1;
		if(s.isgap) continue;
		s.genPdbRecords(&records, ++posi, nid);
	}
	for (auto & r : records) {
		os << r.toString() << std::endl;
	}
}
void NSPproteinrep::readbackbonefrompdb(const std::string &filename, std::vector<std::vector<BackBoneSite>> &chains){
	PdbReader reader;
	reader.readpdb(filename);
	std::string chainids=reader.chainids();
	std::vector<BackBoneSite>  sites;
	for(int i=0;i<chainids.size();++i){
		char chainid=chainids[i];
		std::vector<std::string> seq=reader.getaminoacidsequence(chainid);
		for(int resinumber=0;resinumber<seq.size();++resinumber){
			typename PdbReader::ResKeyType reskey=reader.mappdbkeyint()->pdbResKey(resinumber,i);
			std::vector<PdbRecord> &records=reader.records().at(chainid).at(reskey);
			BackBoneSite bs;
			bs.chainid=chainid;
			bs.resid=records[0].residueid;
			bs.resname=records[0].residuename;
			bs.resseq=resinumber;
			bs.sscode='U';
			std::vector<NSPgeometry::XYZ> crd(4);
			int nbcatoms=0;
			for(auto &r:records) {
				if(r.atomname=="N") {
					crd[0]=NSPgeometry::XYZ(r.x,r.y,r.z);
					++nbcatoms;
				} else if(r.atomname=="CA"){
					crd[1]=NSPgeometry::XYZ(r.x,r.y,r.z);
					++nbcatoms;
				} else if(r.atomname=="C"){
					crd[2]=NSPgeometry::XYZ(r.x,r.y,r.z);
					++nbcatoms;
				} else if(r.atomname=="O"){
					crd[3]=NSPgeometry::XYZ(r.x,r.y,r.z);
					++nbcatoms;
				}
			}
			assert(nbcatoms>=4);
			bs.changecrd(crd);
			sites.push_back(bs);
		}
	}
	chains.assign(1,std::vector<BackBoneSite>());
	int chainnumber=0;
	for(auto it=sites.begin();it !=sites.end()-1;++it){
		chains[chainnumber].push_back(*it);
		const double *d0 = it->data();
		NSPgeometry::XYZ c0(d0[BackBoneSite::CCRD], d0[BackBoneSite::CCRD + 1],
				d0[BackBoneSite::CCRD + 2]);
		const double *d1 = (it+1)->data();
		NSPgeometry::XYZ n1(d1[BackBoneSite::NCRD], d1[BackBoneSite::NCRD + 1],
				d1[BackBoneSite::NCRD + 2]);
		double bond = NSPgeometry::distance(c0, n1);
		if (bond > 2.0){
			chains.push_back(std::vector<BackBoneSite>());
			++chainnumber;
		}
	}
	chains[chainnumber].push_back(sites.back());
}
std::vector<BackBoneSite> NSPproteinrep::generaterandombackbone(int length,
		const std::vector<std::pair<int,int>> &helixregions,
		const std::vector<std::pair<int,int>> &strandregions){
	auto rng=NSPdstl::RandomEngine<>::getinstance().realrng(0,1);
	std::vector<BackBoneSite> chain(length);
	for(int i=0;i<chain.size();++i){
		bool helix=false;
		bool strand=false;
		for(auto &r:helixregions){
			if(i>=r.first &&i<r.second) {
				helix=true;break;
			}
		}
		for(auto &r:strandregions){
			if(i>=r.first &&i<r.second) {
				strand=true;break;
			}
		}
		double phi,psi;
		if(helix) NSPpdbstatistics::PhiPsiDistr::helixdistr().randomphipsi(rng,&phi,&psi);
		else if(strand)NSPpdbstatistics::PhiPsiDistr::stranddistr().randomphipsi(rng,&phi,&psi);
		else NSPpdbstatistics::PhiPsiDistr::mixcoildistr().randomphipsi(rng,&phi,&psi);
		BackBoneSite *bsp=nullptr;
		if(i>0) bsp=&chain[i-1];
		genbackbonesite(bsp,false,phi,psi,&chain[i]);
		if(i==0) {
			chain[0].chainid='A';
			chain[0].resseq=0;
		}
		if(helix) chain[i].sscode='H';
		else chain[i].sscode='C';
	}
	return chain;
}
std::vector<int> NSPproteinrep::findcissites(const std::vector<BackBoneSite> &chain){
	std::vector<int> res;
	for (int i=0;i<chain.size()-1;++i) {
		BackBoneSite bs=chain[i];
		double omiga=bs.omiga(chain[i+1]);
		if(omiga>-90.0 &&omiga<90.0) res.push_back(i+1);
	}
	return res;
}
char BackBoneSite::sscodechar() const {
	char c = sscode;
	if (c != 'H' && c != 'E')
		return 'C';
	return c;
}
NSPgeometry::LocalFrame BackBoneSite::localframe() {
	return NSPgeometry::make_localframe(getcrd(CACRD), getcrd(NCRD),
			getcrd(CCRD));
}
BBatomdis2::BBatomdis2(BackBoneSite &s1, BackBoneSite &s2) {
	std::vector<NSPgeometry::XYZ> crd1;
	std::vector<NSPgeometry::XYZ> crd2;
	crd1.push_back(s1.hcrd());
	s1.getcrd(crd1);
	crd2.push_back(s2.hcrd());
	s2.getcrd(crd2);
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 5; ++j) {
			dist2_.push_back(NSPgeometry::distance2(crd1[i], crd2[j]));
		}
	}
}

NSPgeometry::QuaternionCrd NSPproteinrep::relativeorientation(BackBoneSite &s1,
		BackBoneSite &s2) {
	NSPgeometry::QuaternionCrd Q1(s1.localframe()); //global s1 to local
	NSPgeometry::QuaternionCrd Q2(s2.localframe()); //global s2 to local
	return (Q1.invert()) * Q2; //global s2 to global s1
}
NSPgeometry::RelativePosition NSPproteinrep::relativeposition(BackBoneSite &s1,
		BackBoneSite &s2) {
	NSPgeometry::LocalFrame lf1 = s1.localframe();
	NSPgeometry::XYZ ca2 = s2.getcrd(BackBoneSite::CACRD);
	ca2 = lf1.global2localcrd(ca2);
	NSPgeometry::QuaternionCrd Q = relativeorientation(s1, s2);
	NSPgeometry::RelativePosition rp;
	rp.location = ca2;
	rp.orientation = Q;
	return rp;
}
int NSPproteinrep::nextresidue_atomsclashed(const BackBoneSite &s1,
		const BackBoneSite &s2) {
	double rmin2_1=8.12;
	double rmin2_2=9.0;
//	double distoc2 = NSPgeometry::distance2(s1.ocrd(), s2.ccrd());
//	double rmin2=distoc2;
	double distoo2 = NSPgeometry::distance2(s1.ocrd(), s2.ocrd());
	double rmin2=distoo2;
	double distnca2 = NSPgeometry::distance2(s1.ncrd(), s2.cacrd());
	rmin2=distnca2<rmin2?distnca2:rmin2;
	double distnc2 = NSPgeometry::distance2(s1.ncrd(), s2.ccrd());
	rmin2=distnc2<rmin2?distnc2:rmin2;
	double distno2 = NSPgeometry::distance2(s1.ncrd(), s2.ocrd());
	if (distno2 < 5.29){
//		std::cout<< "distno2 " << distno2<<std::endl;
		return 2;
	}
	double distcac2 = NSPgeometry::distance2(s1.cacrd(), s2.ccrd());
	rmin2=distcac2<rmin2?distcac2:rmin2;
	double distcao2 = NSPgeometry::distance2(s1.cacrd(), s2.ocrd());
	rmin2=distcao2<rmin2?distcao2:rmin2;
	double distco2 = NSPgeometry::distance2(s1.ccrd(), s2.ocrd());
	rmin2=distco2<rmin2?distco2:rmin2;
	if(rmin2<rmin2_1) return 3;
	if(rmin2 <rmin2_2) return 1;
	return 0;
}
bool NSPproteinrep::atomsclashed(const BackBoneSite &s1,
		const BackBoneSite &s2) {
	for (int i = 0; i < 4; i++) {
		int icrd = BackBoneSite::NCRD + 3 * i;
		NSPgeometry::XYZ ri = s1.getcrd(icrd);
		for (int j = 0; j < 4; j++) {
			int jcrd = BackBoneSite::NCRD + 3 * j;
			NSPgeometry::XYZ rj = s2.getcrd(jcrd);
			double dis2 = NSPgeometry::distance2(ri, rj);
			if (dis2 < 6.25)
				return true;
			if (dis2 < 9
					&& (!((icrd == BackBoneSite::NCRD
							&& jcrd == BackBoneSite::OCRD)
							|| (icrd == BackBoneSite::OCRD
									&& jcrd == BackBoneSite::NCRD))))
				return true;
		} //j
	} //i
	return false;
}

void NSPproteinrep::readbackbonesites(const std::string & filename,
		std::vector<BackBoneSite> &sites) {
	std::ifstream is;
	is.open(filename.c_str());
	if (is.fail()) {
		std::cout << "Cannot open sites data file: " << filename << std::endl;
		std::cout << "Program will exit." << std::endl;
		exit(0);
	}
	char buffer[120];
	int nsection = 0;
	while (is.getline(buffer, 120)) {
		std::vector<std::string> bslines;
		bslines.push_back(std::string(buffer));
		for (int i = 0; i < 4; i++) {
			is.getline(buffer, 120);
			bslines.push_back(std::string(buffer));
		}
		auto it = bslines.begin();
		try {
			sites.push_back(BackBoneSite::read(it));
		} catch (std::exception &err) {
			std::cout << "Error encountered reading backbone sites "
					<< std::endl;
			exit(1);
		}
		if (++nsection == 10000) {
			nsection = 0;
			std::cout << "number of sites read: " << sites.size() << std::endl;
		}
	}
	is.close();
}

bool NSPproteinrep::readbackbonesites(std::istream &is, int number,
		std::vector<BackBoneSite> &sites) {
	char buffer[120];
	for (int ns = 0; ns < number; ++ns) {
		std::vector<std::string> bslines;
		for (int i = 0; i < 5; i++) {
			is.getline(buffer, 120);
			if (!is.good())
				return false;
			bslines.push_back(std::string(buffer));
		}
		auto it = bslines.begin();
		try {
			sites.push_back(BackBoneSite::read(it));
		} catch (std::exception &err) {
			return false;
		}
	}
	return true;
}
void NSPproteinrep::genbackbonesite(BackBoneSite *psite, bool cispep,
		double phi, double psi, BackBoneSite *newsite) {
	if (psite) {
		newsite->pdbid = psite->pdbid;
		newsite->chainid = psite->chainid;
		newsite->resid = psite->resid + 1;
		newsite->resseq = psite->resseq + 1;
		if (cispep)
			psite->data_[BackBoneSite::OMIGA] = 0.0;
	} else {
		newsite->resid = 1;
		newsite->resseq = 1;
	}
	while (phi > 180.0)
		phi -= 360.0;
	while (phi < -180.0)
		phi += 360.0;
	while (psi > 180.0)
		psi -= 360.0;
	while (psi < -180.0)
		psi += 360.0;
	newsite->data_[BackBoneSite::PHI] = phi;
	newsite->data_[BackBoneSite::PSI] = psi;
	newsite->data_[BackBoneSite::OMIGA] = 180.0;
	std::vector<NSPgeometry::XYZ> crd;
	IdealGeometries & igdat = IdealGeometries::getGlobalInstance();
	double degree = 3.14159265 / 180.0;
	if (psite) {
		crd.push_back(
				NSPgeometry::InternaltoXYZ(psite->ccrd(), psite->cacrd(),
						psite->ncrd(), igdat.idealLength("C", "N"),
						igdat.idealAngle("CA", "C", "N"),
						psite->psi() * degree));
		crd.push_back(
				NSPgeometry::InternaltoXYZ(crd[0], psite->ccrd(),
						psite->cacrd(), igdat.idealLength("N", "CA"),
						igdat.idealAngle("C", "N", "CA"),
						psite->omiga() * degree));
		crd.push_back(
				NSPgeometry::InternaltoXYZ(crd[1], crd[0], psite->ccrd(),
						igdat.idealLength("CA", "C"),
						igdat.idealAngle("N", "CA", "C"), phi * degree));
	} else {
		crd.push_back(NSPgeometry::XYZ());
		crd.push_back(
				NSPgeometry::InternaltoXYZ(crd[0],
						igdat.idealLength("N", "CA")));
		crd.push_back(
				NSPgeometry::InternaltoXYZ(crd[1], crd[0],
						igdat.idealLength("CA", "C"),
						igdat.idealAngle("N", "CA", "C")));
	}
	crd.push_back(
			NSPgeometry::InternaltoXYZ(crd[2], crd[1], crd[0],
					igdat.idealLength("C", "O"),
					igdat.idealAngle("CA", "C", "O"), (psi + 180.0) * degree));
	newsite->changecrd(crd);
}
void NSPproteinrep::genprevbackbonesite(BackBoneSite *nsite,double omiga,double psi, double phi,
		BackBoneSite *newsite){
		newsite->pdbid = nsite->pdbid;
		newsite->chainid = nsite->chainid;
		newsite->resid = nsite->resid - 1;
		newsite->resseq = nsite->resseq - 1;
		while (phi > 180.0)
			phi -= 360.0;
		while (phi < -180.0)
			phi += 360.0;
		while (psi > 180.0)
			psi -= 360.0;
		while (psi < -180.0)
			psi += 360.0;
		newsite->data_[BackBoneSite::PHI] = phi;
		newsite->data_[BackBoneSite::PSI] = psi;
		newsite->data_[BackBoneSite::OMIGA] = omiga;
		std::vector<NSPgeometry::XYZ> crd;
		IdealGeometries & igdat = IdealGeometries::getGlobalInstance();
		double degree = 3.14159265 / 180.0;
		crd.resize(4);
		crd[2]=NSPgeometry::InternaltoXYZ(nsite->ncrd(), nsite->cacrd(),
						nsite->ccrd(), igdat.idealLength("C", "N"),
						igdat.idealAngle("C", "N", "CA"),
						nsite->phi() * degree);
		crd[3]=NSPgeometry::InternaltoXYZ(crd[2], nsite->ncrd(),
						nsite->cacrd(), igdat.idealLength("C", "O"),
						igdat.idealAngle("O", "C", "N"),
						(180.0-omiga) * degree);
		crd[1]=NSPgeometry::InternaltoXYZ(crd[2], nsite->ncrd(), nsite->cacrd(),
						igdat.idealLength("CA", "C"),
						igdat.idealAngle("CA", "C", "N"), omiga * degree);
		crd[0]=	NSPgeometry::InternaltoXYZ(crd[1], crd[2], nsite->ncrd(),
					igdat.idealLength("N", "CA"),
					igdat.idealAngle("N", "CA", "C"), psi * degree);
		newsite->changecrd(crd);
}
std::vector<BackBoneSite> NSPproteinrep::loops2gaps(const std::vector<BackBoneSite> & chain,
		int extension){
	std::vector<BackBoneSite> newchain(chain.size());
	std::copy(chain.begin(),chain.end(),newchain.begin());
	for(auto &s:newchain) s.isgap=true;
	for(int i=0;i<newchain.size();++i){
		char sscode=newchain[i].sscode;
		if(sscode=='H'||sscode =='E'||sscode=='m'||sscode=='d') {
			newchain[i].isgap=false;
			for(int offset=0;offset<extension;++offset){
				if(i>offset) newchain[i-offset-1].isgap=false;
				if(i<newchain.size()-1-offset) newchain[i+offset+1].isgap=false;
			}
		}
	}
	return newchain;
}
std::vector<BackBoneSite> NSPproteinrep::insertgaps(const std::vector<BackBoneSite> & chain,
		const std::map<int,int> & gapposi_lens){
	std::vector<BackBoneSite> newchain;
	int resseq=0;
	for(int i=0;i<chain.size();++i){
		if(gapposi_lens.find(i) != gapposi_lens.end()){
			int len=gapposi_lens.at(i);
			for(int l=0;l<len; ++i){
				newchain.push_back(BackBoneSite());
				newchain.back().isgap=true;
				newchain.back().resseq=resseq++;
			}
		}
		newchain.push_back(chain[i]);
		newchain.back().resseq=resseq++;
	}
	return newchain;
}
NSPdstl::AlignedPositions NSPproteinrep::alignaspointsets(
		const std::vector<BackBoneSite> &seta, const std::vector<BackBoneSite> &setb,
		double rcut){
	std::vector<NSPgeometry::XYZ> crda;
	std::vector<NSPgeometry::XYZ> crdb;
	for(auto &s:seta) crda.push_back(s.cacrd());
	for(auto &s:setb) crdb.push_back(s.cacrd());
	return NSPgeometry::alignpointset(crda,crdb,rcut);
}
int NSPproteinrep::sselements(const std::vector<BackBoneSite> &chain, std::vector<int> *elmtidseq,
		std::vector<std::pair<int,int>> *elmtposlens){
	elmtidseq->clear();
	elmtidseq->resize(chain.size(),-1);
	elmtposlens->clear();
	bool inss=false;
	const BackBoneSite *sprev=&(chain[0]);
	int startposi=-1;
	int len;
	if(sprev->sscodechar()!='C') {
		inss=true;
		startposi=0;
		len=1;
	}
	for(int i=1;i<chain.size();++i){
		const BackBoneSite &s=chain[i];
		bool consecutive=((s.ncrd()-sprev->ccrd()).squarednorm()<3.24);
		if(!consecutive) {
			if(inss){
				elmtposlens->push_back(std::make_pair(startposi,len));
				inss=false;
			}
		}
		if(!inss){
			if(s.sscodechar() !='C') {
				startposi=i;
				inss=true;
				len=1;
			}
		} else {
			if(s.sscodechar()!='C'){
				if(s.sscodechar()==sprev->sscodechar()){
					len++;
				} else {
					elmtposlens->push_back(std::make_pair(startposi,len));
					startposi=i;
					len=1;
				}
			} else {
				elmtposlens->push_back(std::make_pair(startposi,len));
				inss=false;
			}
		}
		sprev=&s;
	}
	if(inss) elmtposlens->push_back(std::make_pair(startposi,len));
	int elmtid=0;
	for(auto &elmt:*elmtposlens){
		for(int i=0;i<elmt.second;++i){
			elmtidseq->at(elmt.first+i)=elmtid;
		}
		++elmtid;
	}
	return elmtposlens->size();
}
std::vector<double> NSPproteinrep::extractcrd(const std::vector<NSPproteinrep::BackBoneSite> &chain){
	std::vector<double> crd;
	for(auto iter=chain.begin(); iter !=chain.end();++iter){
		for(int d=0;d<3;++d) crd.push_back(iter->ncrd()[d]);
		for(int d=0;d<3;++d) crd.push_back(iter->cacrd()[d]);
		for(int d=0;d<3;++d) crd.push_back(iter->ccrd()[d]);
		for(int d=0;d<3;++d) crd.push_back(iter->ocrd()[d]);
	}
	return crd;
}
void NSPproteinrep::assigncrd(const std::vector<double> &crd, std::vector<BackBoneSite> &chain){
	assert(crd.size()==12*chain.size());
	int idx=0;
	for(auto iter=chain.begin(); iter !=chain.end();++iter){
		std::vector<double> crdi(12);
		std::copy(crd.begin()+12*idx, crd.begin()+12*(idx+1), crdi.begin());
		iter->changecrd(crdi);
		++idx;
	}
	for(auto iter=chain.begin(); iter !=chain.end();++iter){
		if(iter!=chain.begin())
			iter->phi(*(iter-1));
		if((iter+1) != chain.end()){
			iter->psi(*(iter+1));
			iter->omiga(*(iter+1));
		}
	}
}
/*
 void pdbio::make_backbonesites_chain(ChainData & chain,
 ChainAttrib<BackBoneSite::backbonesite_attrib_type> &sites){
 PdbStructure::Chain *c=chain.pdb->getchain(chain.chainid);
 if(!c->has_attrib(StrideReader::attrib_name)) return;
 StrideReader::stride_attrib_type *strides=
 c->get_attrib<StrideReader::stride_attrib_type>(StrideReader::attrib_name);

 int ir=0;
 for (auto & res : c->residues) {
 posidx_t N0=locate_atom(*chain.pdb,
 res.first,
 res.second,
 "N");
 posidx_t CA0=locate_atom(*chain.pdb,
 res.first,
 res.second,
 "CA");
 posidx_t C0=locate_atom(*chain.pdb,
 res.first,
 res.second,
 "C");
 posidx_t O0=locate_atom(*chain.pdb,
 res.first,
 res.second,
 "O");
 //		assert(N0>=0 && CA0>=0 && C0>=0 && O0 >=0);
 if(!(N0>=0 && CA0>=0 && C0>=0 && O0 >=0)) {
 std::cout <<chain.pdbid << chain.chainid
 <<(*chain.pdb)[res.first].residueid <<std::endl;
 abort();
 }
 double omiga=0.0;
 if(ir+1< c->residues.size()) {
 posidx_t N1=locate_atom(*chain.pdb,c->residues[ir+1].first,
 c->residues[ir+1].second,"N");
 posidx_t CA1=locate_atom(*chain.pdb,c->residues[ir+1].first,
 c->residues[ir+1].second,"CA");
 assert(N1>=0 and CA1 >=0);
 double bond=NSPgeometry::distance((*chain.pdb)[C0],(*chain.pdb)[N1]);
 if(bond < 2.0) omiga=
 NSPgeometry::torsion((*chain.pdb)[CA0],
 (*chain.pdb)[C0],
 (*chain.pdb)[N1],
 (*chain.pdb)[CA1]);
 };
 BackBoneSite site;
 site.pdbid=chain.pdbid;
 site.chainid=chain.chainid;
 site.resid=(*chain.pdb)[N0].residueid;
 site.resname=(*chain.pdb)[N0].residuename;
 site.resseq=ir;
 //		assert (site.resid == strides->at(ir).resid);
 if(!site.resid == strides->at(ir).resid) {
 //		std::cout <<strides->at(ir).resid <<"\n";
 //		if (ir >= strides->size()){
 std::cout <<ir <<" "<<strides->size() <<"\n";
 std::cout <<strides->at(ir-1).resid <<"\n";
 std::cout <<chain.pdbid << chain.chainid
 <<(*chain.pdb)[res.first].residueid <<std::endl;
 abort();
 }
 site.sscode=strides->at(ir).sscode;
 site.data()[BackBoneSite::PHI]=strides->at(ir).phi;
 site.data()[BackBoneSite::PSI]=strides->at(ir).psi;
 site.data()[BackBoneSite::OMIGA]=omiga*180.0/3.14159265;
 site.data()[BackBoneSite::SASA]=strides->at(ir).sasa;
 site.data()[BackBoneSite::NCRD]=(*chain.pdb)[N0].x;
 site.data()[BackBoneSite::NCRD+1]=(*chain.pdb)[N0].y;
 site.data()[BackBoneSite::NCRD+2]=(*chain.pdb)[N0].z;
 site.data()[BackBoneSite::CACRD]=(*chain.pdb)[CA0].x;
 site.data()[BackBoneSite::CACRD+1]=(*chain.pdb)[CA0].y;
 site.data()[BackBoneSite::CACRD+2]=(*chain.pdb)[CA0].z;
 site.data()[BackBoneSite::CCRD]=(*chain.pdb)[C0].x;
 site.data()[BackBoneSite::CCRD+1]=(*chain.pdb)[C0].y;
 site.data()[BackBoneSite::CCRD+2]=(*chain.pdb)[C0].z;
 site.data()[BackBoneSite::OCRD]=(*chain.pdb)[O0].x;
 site.data()[BackBoneSite::OCRD+1]=(*chain.pdb)[O0].y;
 site.data()[BackBoneSite::OCRD+2]=(*chain.pdb)[O0].z;
 sites.attrib()->push_back(site);
 //		std::cout<<site.toString();
 ir++;
 }
 sites.attach(c,BackBoneSite::attrib_name);
 //	sites.chain()=c;
 //	c->add_attrib(BackBoneSite::attrib_name,sites.attrib());
 }

 BackBoneHB pdbio::hbproperty(BackBoneSite &s1, BackBoneSite &s2){
 NSPgeometry::XYZ o1=s1.getcrd(BackBoneSite::OCRD);
 NSPgeometry::XYZ n2=s2.getcrd(BackBoneSite::NCRD);
 double dis2_12=NSPgeometry::distance2(o1,n2);
 NSPgeometry::XYZ n1=s1.getcrd(BackBoneSite::NCRD);
 NSPgeometry::XYZ o2=s2.getcrd(BackBoneSite::OCRD);
 double dis2_21=NSPgeometry::distance2(n1,o2);
 BackBoneSite *don=&s1;
 BackBoneSite *acc=&s2;
 if(dis2_12 < dis2_21){
 don=&s2;acc=&s1;
 }
 NSPgeometry::XYZ rx=acc->getcrd(BackBoneSite::CCRD);
 NSPgeometry::XYZ ra=acc->getcrd(BackBoneSite::OCRD);
 NSPgeometry::XYZ rh=don->hcrd();
 NSPgeometry::XYZ rd=don->getcrd(BackBoneSite::NCRD);
 return BackBoneHB(acc->sscodechar(),don->sscodechar(),rx,ra,rh,rd);
 }
 */

