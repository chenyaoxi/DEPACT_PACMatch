/*
 * ProteinRep.h
 *
 *  Created on: 2017年10月24日
 *      Author: notxp
 */

#ifndef DESIGNSEQ_PROTEINREP_H_
#define DESIGNSEQ_PROTEINREP_H_
#include <string>
#include "stdio.h"
#include <vector>
#include <map>
#include <set>
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "geometry/quatfit.h"
#include "designseq/StringTool.h"
#include "designseq/AtomLib.h"
#include "designseq/Rotamer.h"
#include "backbone/backbonesite.h"

namespace NSPdesignseq {
using namespace std;
using namespace NSPgeometry;
using namespace NSPproteinrep;

class Atom{
public:
	string resType;
	string name;
	string type;
	XYZ coord;

	Atom();
	Atom(string line);
	Atom(string name, XYZ coord);
	Atom& operator=(const Atom& other);
	bool isBackboneAtom() const;
	void guessAtomType();
	void setAtomType(string type);
	void setResType(string resType);
	bool isIon() const;
	string& getName();
	string& getType();
	XYZ& getCoord();
	string& getResType();
	void setCoord(const XYZ& coord);
	float distance(const Atom& other) const;
	string nameString()	 const;
	virtual ~Atom();
};

class Residue{
private:
	vector<Atom*> atomList;
	vector<Atom*> backboneAtoms;
	vector<Atom*> sidechainAtoms;
	map<string, Atom*> atomMap;
public:
	string resID;
	int resSeqID;
	string triName;
	char chainID;
	int atomNum;
	bool hasLocalFrame;
	LocalFrame coordSys;
	char altLoc;
	Residue();
	Residue(string resID, char chainID, string triName);
	Residue(NSPproteinrep::BackBoneSite* bbSite);
	void addAtom(Atom* a);
	void setResSeqID(int id);
	void setAltLoc(char c) {this->altLoc = c;}
	bool hasThreeCoreAtoms() const;
	void updateCoordSystem();
	bool hasAtom(const string& atomName) const;
	Atom* getAtom(const string& atomName);
	vector<Atom*>* getAtomList();
	vector<Atom*>* getBackboneAtoms();
	vector<Atom*>* getSidechainAtoms();
	bool sidechainComplete(AtomLib* atomLib) const;
	XYZ getCbCoord();
	LocalFrame& getCoordSystem();
	char getChainID() const;
	int getResSeqID() const;
	string getResID() const;
	string getType() const;
	int buildRotamer(Rotamer* rot);
	NSPproteinrep::BackBoneSite getBackboneSite();



	Rotamer natRotamer(AtomLib* atLib);

	int printPDBFormat(ofstream& out, int startAtomID) const;
	virtual ~Residue();
};

class PolarAtom{
private:
	string uniqueName;
	XYZ support;
	XYZ core;
	bool isDonor;
	bool isAcceptor;
	float vdwRadius;
public:
	PolarAtom(Residue* res, string atomName);
	PolarAtom(Residue* res, string atomName, Atom* preC);
	PolarAtom(XYZ core, XYZ support, AtomProperty* ap){
		this->core = core;
		this->support = support;
		this->uniqueName = ap->atomUniqueName;
		this->isDonor = ap->isHDonor;
		this->isAcceptor = ap->isHAcceptor;
		this->vdwRadius = ap->vdwRadius;
	}
	XYZ& getCore() {return this->core;}
	XYZ& getSupport() {return this->support;}
	bool isDonerAtom() {return this->isDonor;}
	bool isAcceptorAtom() {return this->isAcceptor;}
	bool neighborTo(PolarAtom* other) {return this->core.distance(other->core) < (this->vdwRadius + other->vdwRadius + 2.0);}
	string getName() {return this->uniqueName;}
	float getVdwRadius() {return this->vdwRadius;}

	virtual ~PolarAtom();
};

class ProteinChain{
private:
	string pdbID;
	char chainID;
	int chainLen;
	vector<Residue*> resList;
	map<string,Residue*> resMap;

public:
	ProteinChain();
	ProteinChain(string pdbID, char chainID);
	ProteinChain(char chainID);
	ProteinChain(vector<NSPproteinrep::BackBoneSite>* bbSites);

	void setPDBID(string pdbID);
	void setChainID(char c);
	string getPDBID() const;
	char getChainID() const;
	int getChainLength() const;
	vector<Residue*>& getResList();
	Residue* getResidue(const string& resID);
	void addResidue(Residue* res);
	string getSequence() const;
	int printPDBFormat(ofstream& out, int startAtomID) const;
	int printPDBFormatNoHydrogen(ofstream& out, int startAtomID) const;
	virtual ~ProteinChain();
};

class PDB{
private:
	string pdbID;
	vector<ProteinChain*> chains;
	vector<Residue*> residues;
public:
	PDB();
	PDB(const string& pdbFile, const string& pdbID);
	void init(std::istream &is);
	PDB(std::istream &is, const string &pdbID){
		this->pdbID=pdbID;
		init(is);
	}
	PDB& operator=(const PDB& other);
	vector<ProteinChain*>& getChains();
	ProteinChain* getFirstChain();
	ProteinChain* getChain(char c);
	string getFirstSeq();
	vector<Residue*>& getResList();
	vector<Residue*> getValidResList(){
		vector<Residue*> list;
		for(size_t i=0;i<residues.size();i++){
			Residue* res = residues.at(i);
			if(res->hasThreeCoreAtoms())
				list.push_back(res);
		}
		return list;
	}
	void printPDBFormat(ofstream& out) const;
	void printPDBFormatNoHydrogen(ofstream& out) const;
	string getPDBID()
	{
		return this->pdbID;
	}
	virtual ~PDB();
};

class Phipsi{
public:
	float phi;
	float psi;
	Phipsi() {this->phi = 0; this->psi = 0;}
	Phipsi(float phi, float psi){
		// TODO Auto-generated constructor stub
		this->phi = phi;
		this->psi = psi;
		if(phi == 360)
			this->phi = -60;
		if(psi == 360)
			this->psi = 130;

		if(this->phi >= 180)
			this->phi = this->phi - 360;
		if(this->phi < -180)
			this->phi = this->phi + 360;
		if(this->psi >= 180)
			this->psi = this->psi - 360;
		if(this->psi < -180)
			this->psi = this->psi + 360;

		if(this->phi > 180 || this->phi < -180 || this->psi > 180 || this->psi < -180)
		{
			cerr << "invalid phi psi: " << phi << " " << psi << endl;
			exit(1);
		}
	}
	float distance(const Phipsi& other) const;
	char regionAB() const;
	virtual ~Phipsi();
};

class PhipsiLib{
private:
	int pointNum;
	vector<Phipsi*> ppList;
	int indexTable[36][36][20];
	void creatIndexTable();
public:
	PhipsiLib();
	Phipsi* indexToPhipsi(int id) const;
	int phipsiToIndex(const Phipsi* pp) const;
	int findNearestPointsWithoutIndex(Phipsi* pp);
	virtual ~PhipsiLib();
};

class ResPairOrientation{
private:

	vector<XYZ> points;

public:
	ResPairOrientation();
	ResPairOrientation(LocalFrame& csA, LocalFrame& csB);
	ResPairOrientation(LocalFrame& csB);
	ResPairOrientation(NSPproteinrep::BackBoneSite& resA, NSPproteinrep::BackBoneSite& resB);
	ResPairOrientation(string& s);

	void addAtoms(LocalFrame& csB);
	double rmsd(const ResPairOrientation& other) const;
	XYZ getTern(int index);
	string toString() const;
};

class SaiPair{
public:
	float saiA;
	float saiB;

	SaiPair(float x, float y){
		this->saiA = x;
		this->saiB = y;
	}

	float distanceSquare(SaiPair* other){
		float dx = this->saiA - other->saiA;
		float dy = this->saiB - other->saiB;
		return dx*dx+dy*dy;
	}
};

class BackboneSitesPair{
public:
	NSPproteinrep::BackBoneSite* siteA;
	NSPproteinrep::BackBoneSite* siteB;
	ResPairOrientation ori;
	int seqSep;

	BackboneSitesPair(BackBoneSite* siteA, BackBoneSite* siteB);


};

LocalFrame getBackboneSiteLocalFrame(const NSPproteinrep::BackBoneSite& bbSite);
LocalFrame buildLocalFrame(LocalFrame& csNewGlobal, LocalFrame& csOld);

} /* namespace NSPdesignseq*/




#endif /* DESIGNSEQ_PROTEINREP_H_ */
