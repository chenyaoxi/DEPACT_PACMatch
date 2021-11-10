/*
 * ProteinRep.cpp
 *
 *  Created on: 2017年10月24日
 *      Author: notxp
 */

#include "designseq/ProteinRep.h"

namespace NSPdesignseq {

Atom::Atom(){
	this->name = "X";
	this->type = "X";
	this->resType = "";
	this->coord = XYZ();
}

Atom::Atom(string line){
	string s = line.substr(12,4);
	this->name = trimString(s);
	this->type = "";
	float x = atof(line.substr(30,8).c_str());
	float y = atof(line.substr(38,8).c_str());
	float z = atof(line.substr(46,8).c_str());
	this->coord = XYZ(x,y,z);
	guessAtomType();
	s = line.substr(17,3);
	this->resType = trimString(s);
}

Atom::Atom(string name, XYZ coord){
	this->name = name;
	this->coord = coord;
	this->type = "X";
	this->resType = "UNK";
}


Atom& Atom::operator=(const Atom& other)
{
	if(this == &other)
		return *this;
	this->name = other.name;
	this->type = other.type;
	this->coord = other.coord;
	this->resType = other.resType;
	return *this;
}

bool Atom::isBackboneAtom() const
{
	return name=="N" || name=="CA" || name=="C" || name=="O" || name=="OXT" || name=="OT1" || name=="OT2";
}

void Atom::guessAtomType()
{
	char c = name.at(0);
	if(name.length() == 1)
		this->type = name;
	else if(name == "FE" || name == "ZN" || name=="MG" || name=="MN" || name=="CL" || name == "SE" || name=="CU" || name=="NA")
		this->type = name;
	else if(c>='A' && c <= 'Z')
	{
	    this->type = string(1,c);

	}
	else
		this->type = string(1,name.at(1));
}

bool Atom::isIon() const
{
	if(type == "FE" || type == "ZN" || type == "MG" || type == "MN" || type == "CA" || type == "CU" || type == "NI")
		return true;
	return false;
}

string& Atom::getName()
{
	return this->name;
}

string& Atom::getType()
{
	return this->type;
}

string& Atom::getResType()
{
	return this->resType;
}



XYZ& Atom::getCoord()
{
	return this->coord;
}

void Atom::setCoord(const XYZ& coord)
{
	this->coord = coord;
}

float Atom::distance(const Atom& other) const
{
	return this->coord.distance(other.coord);
}

void Atom::setAtomType(string type)
{
	this->type = type;
}

void Atom::setResType(string resType)
{
	this->resType = resType;
}

string Atom::nameString() const
{
	char s[10];
	if(this->type.length() > 1)
		sprintf(s,"%-5s",this->name.c_str());
	else
		sprintf(s," %-4s",this->name.c_str());
	string ss = string(s);
	return ss;
}

Atom::~Atom(){

}

Residue::Residue() {
	this->resID = "1";
	this->resSeqID = 0;
	this->atomNum = 0;
	this->triName = "UNK";
	this->chainID = '-';
	this->hasLocalFrame = false;
	this->coordSys;
	this->altLoc = ' ';
	// TODO Auto-generated constructor stub

}

Residue::Residue(string resID, char chainID, string triName)
{
	this->resID = resID;
	this->chainID = chainID;
	this->triName = triName;
	this->resSeqID = 0;
	this->atomNum = 0;
	this->hasLocalFrame = false;
	this->altLoc = ' ';
}

void Residue::addAtom(Atom* a)
{
	this->atomList.push_back(a);
	if(a->isBackboneAtom())
		this->backboneAtoms.push_back(a);
	else
		this->sidechainAtoms.push_back(a);
	this->atomMap[a->getName()] = a;
	this->atomNum ++;
}

void Residue::setResSeqID(int id)
{
    this->resSeqID = id;
}

bool Residue::hasThreeCoreAtoms() const
{
	if(atomMap.find("N") == atomMap.end())
		return false;
	if(atomMap.find("CA") == atomMap.end())
		return false;
	if(atomMap.find("C") == atomMap.end())
		return false;
	return true;
}

void Residue::updateCoordSystem()
{
	if(!this->hasThreeCoreAtoms())
		return;
	Atom* N = this->atomMap.at("N");
	Atom* C = this->atomMap.at("C");
	Atom* CA = this->atomMap.at("CA");

	XYZ n = N->getCoord();
	XYZ c = C->getCoord();
	XYZ ca = CA->getCoord();

	XYZ can = ~(N->getCoord() - CA->getCoord());
	XYZ cac = ~(C->getCoord() - CA->getCoord());

	XYZ z = ~(can^cac);
	XYZ x = ~(can+cac);

	XYZ y = ~(z^x);

	this->coordSys = LocalFrame();
	this->coordSys.origin_ = ca;
	this->coordSys.axis_.push_back(x);
	this->coordSys.axis_.push_back(y);
	this->coordSys.axis_.push_back(z);

}

bool Residue::hasAtom(const string& atomName) const
{
    if(atomMap.find(atomName) == atomMap.end())
		return false;
	else
        return true;
}

Atom* Residue::getAtom(const string& atomName)
{
    map<string,Atom*>::const_iterator it = atomMap.find(atomName);
    if(it == atomMap.end())
        return NULL;
    else
        return it->second;
}

vector<Atom*>* Residue::getAtomList()
{
    return &this->atomList;
}

vector<Atom*>* Residue::getBackboneAtoms()
{
    return &this->backboneAtoms;
}

vector<Atom*>* Residue::getSidechainAtoms()
{
    return &this->sidechainAtoms;
}

bool Residue::sidechainComplete(AtomLib* atomLib) const
{
	vector<string>* scAtoms = atomLib->getAminoAcidSidechainAtomNames(this->triName);
	if(scAtoms == NULL) return false;
	vector<string>::iterator it;
	for(it=scAtoms->begin();it<scAtoms->end();it++)
	{
		string s = *it;
		if(!hasAtom(s))
			return false;
	}
	return true;
}

XYZ Residue::getCbCoord()
{
	if(!hasThreeCoreAtoms())
	{
		cerr << "lack backbone atom information" << endl;
		exit(1);
	}
	if(!this->hasLocalFrame)
		updateCoordSystem();
	XYZ localCB = XYZ(-0.942, 0.009, 1.208);
	return this->coordSys.local2globalcrd(localCB);
}

LocalFrame& Residue::getCoordSystem()
{
	if(!hasThreeCoreAtoms())
	{
		cerr << "lack backbone atom information" << endl;
		exit(1);
	}
	if(!this->hasLocalFrame)
		updateCoordSystem();
	return this->coordSys;
}

char Residue::getChainID() const
{
	return this->chainID;
}

int Residue::getResSeqID() const
{
	return this->resSeqID;
}

string Residue::getResID() const
{
	return this->resID;
}

string Residue::getType() const
{
	return this->triName;
}

int Residue::buildRotamer(Rotamer* rot)
{
	for(unsigned int i=0;i<sidechainAtoms.size();i++)
	{
		Atom* a = sidechainAtoms.at(i);
		atomMap.erase(a->getName());
		delete a;
	}
	sidechainAtoms.clear();
	atomList.clear();
	atomNum = 0;

	this->triName = rot->triName;
	for(unsigned int i=0;i<backboneAtoms.size();i++)
	{
		atomList.push_back(backboneAtoms.at(i));
	}

	LocalFrame& cs = getCoordSystem();
	for(unsigned int i=0;i<rot->atomNameList.size();i++){
		Atom* a = new Atom(rot->atomNameList.at(i), cs.local2globalcrd(rot->coordList.at(i)));
		addAtom(a);
	}

	return 1;
}

NSPproteinrep::BackBoneSite Residue::getBackboneSite(){
	NSPproteinrep::BackBoneSite bs;
	bs.pdbid = "unkn";
	bs.chainid = this->chainID;
	bs.resid = boost::lexical_cast<int>(this->resID);
	bs.resname = this->triName;
	bs.resseq = this->resSeqID;
	bs.sscode = 'X';

	Atom* N = this->getAtom("N");
	Atom* CA = this->getAtom("CA");
	Atom* C = this->getAtom("C");
	Atom* O = this->getAtom("O");
	if(N == NULL || CA == NULL || C == NULL || O == NULL)
	{
		cout << "insufficient backbone atoms: \nchainID: " << this->chainID << "\nresID: " <<   this->resID << endl;
	}

	bs.data_[0] = 360.0; //phi
	bs.data_[1] = 360.0; //psi
	bs.data_[2] = 360.0; //omg
	bs.data_[3] = 0.0; //sasa
	bs.data_[4] = N->getCoord()[0];
	bs.data_[5] = N->getCoord()[1];
	bs.data_[6] = N->getCoord()[2];
	bs.data_[7] = CA->getCoord()[0];
	bs.data_[8] = CA->getCoord()[1];
	bs.data_[9] = CA->getCoord()[2];
	bs.data_[10] = C->getCoord()[0];
	bs.data_[11] = C->getCoord()[1];
	bs.data_[12] = C->getCoord()[2];
	bs.data_[13] = O->getCoord()[0];
	bs.data_[14] = O->getCoord()[1];
	bs.data_[15] = O->getCoord()[2];


	return bs;



}

Rotamer Residue::natRotamer(AtomLib* atLib)
{
	if(!sidechainComplete(atLib))
		return Rotamer();
	if(triName == "GLY")
		return Rotamer("GLY-0","GLY");
	LocalFrame& cs = getCoordSystem();
	Rotamer rot(triName+"-0", triName);
	vector<string>* names = atLib->getAminoAcidSidechainAtomNames(triName);
	for(int i=0;i<names->size();i++)
	{
		Atom* a = getAtom(names->at(i));
		rot.addAtom(a->getName(), cs.global2localcrd(a->getCoord()));
	}
	rot.updateLawOfRing();
	return rot;
}



int Residue::printPDBFormat(ofstream& out, int startAtomID) const
{
    char c = this->resID.at(resID.length()-1);
    char s[100];
    vector<Atom*>::const_iterator it;
    int atomID = startAtomID;
    if(c >= '0' && c <= '9')
    {
        for(it=atomList.begin();it!=atomList.end();it++)
        {
            XYZ& coord = (*it)->getCoord();
            sprintf(s,"ATOM%7d %5s%3s %c%4s    %8.3f%8.3f%8.3f",atomID,(*it)->nameString().c_str(),this->triName.c_str(),this->chainID,this->resID.c_str(),coord[0],coord[1],coord[2]);
            out << s << endl;
            atomID++;
        }
    }
    else
    {
        for(it=atomList.begin();it!=atomList.end();it++)
        {
            XYZ& coord = (*it)->getCoord();
            sprintf(s,"ATOM%7d %5s%3s %c%5s   %8.3f%8.3f%8.3f",atomID,(*it)->nameString().c_str(),this->triName.c_str(),this->chainID,this->resID.c_str(),coord[0],coord[1],coord[2]);
            out << s << endl;
            atomID++;
        }
    }
    return atomID;
}

Residue::~Residue() {
	// TODO Auto-generated destructor stub
}

PolarAtom::PolarAtom(Residue* res, string atomName) {
	// TODO Auto-generated constructor stub

	this->uniqueName = res->triName + "-" + atomName;
	this->isDonor = false;
	this->isAcceptor = false;
//	this->HCore = false;

	if(uniqueName == "ASP-OD1" || uniqueName == "ASP-OD2" || uniqueName == "ASN-OD1")
	{
		this->support = res->getAtom("CG")->getCoord();
		this->core = res->getAtom(atomName)->getCoord();
		this->isAcceptor = true;
	}
	else if(uniqueName == "ASN-ND2")
	{
		this->support = res->getAtom("CG")->getCoord();
		this->core = res->getAtom(atomName)->getCoord();
		this->isDonor = true;
	}
	else if(uniqueName == "GLU-OE1" || uniqueName == "GLU-OE2" || uniqueName == "GLN-OE1")
	{
		this->support = res->getAtom("CD")->getCoord();
		this->core = res->getAtom(atomName)->getCoord();
		this->isAcceptor = true;
	}
	else if(uniqueName == "GLN-NE2")
	{
		this->support = res->getAtom("CD")->getCoord();
		this->core = res->getAtom(atomName)->getCoord();
		this->isDonor = true;
	}
	else if(uniqueName == "HIS-ND1")
	{
		XYZ n =  res->getAtom("ND1")->getCoord();
		XYZ sup1 = res->getAtom("CG")->getCoord();
		XYZ sup2 = res->getAtom("CE1")->getCoord();

		this->core = n;
		this->support = sup1 + sup2 -n;

		this->isAcceptor = true;
		this->isDonor = true;
	//	this->HCore = true;

	}
	else if(uniqueName == "HIS-NE2")
	{
		XYZ n = res->getAtom("NE2")->getCoord();

		XYZ sup1 = res->getAtom("CD2")->getCoord();
		XYZ sup2 = res->getAtom("CE1")->getCoord();
		this->core = n;
		this->support = sup1 + sup2 -n;
		this->isAcceptor = true;
		this->isDonor = true;
	//	this->HCore = true;
	}
	else if(uniqueName == "LYS-NZ")
	{
		this->core = res->getAtom(atomName)->getCoord();
		this->support = res->getAtom("CE")->getCoord();
		this->isDonor = true;
	}
	else if(uniqueName == "ARG-NE")
	{
		XYZ n = res->getAtom("NE")->getCoord();
		XYZ sup1 = res->getAtom("CD")->getCoord();
		XYZ sup2 = res->getAtom("CZ")->getCoord();
		this->core = n;
		this->support = sup1 + sup2 -n;
		this->isDonor = true;
	//	this->HCore = true;
	}
	else if(uniqueName == "ARG-NH1" || uniqueName == "ARG-NH2")
	{
		this->core = res->getAtom(atomName)->getCoord();
		this->support = res->getAtom("CZ")->getCoord();
		this->isDonor = true;
	}
	else if(uniqueName == "SER-OG" || uniqueName == "THR-OG1")
	{
		this->core = res->getAtom(atomName)->getCoord();
		this->support = res->getAtom("CB")->getCoord();
		this->isAcceptor = true;
		this->isDonor = true;
	}
	else if(uniqueName == "TRP-NE1")
	{
		XYZ n = res->getAtom("NE1")->getCoord();
		XYZ sup1 = res->getAtom("CD1")->getCoord();
		XYZ sup2 = res->getAtom("CE2")->getCoord();
		this->core = n;
		this->support = sup1 + sup2 -n;
		this->isDonor = true;
	//	this->HCore = true;
	}
	else if(uniqueName == "TYR-OH")
	{
		this->core = res->getAtom(atomName)->getCoord();
		this->support = res->getAtom("CZ")->getCoord();
		this->isAcceptor = true;
		this->isDonor = true;
	}
	else if(atomName == "O" || atomName == "OXT" || atomName == "OT1" || atomName == "OT2")
	{
		this->uniqueName = res->triName + "-O";
		this->core = res->getAtom(atomName)->getCoord();
		this->support = res->getAtom("C")->getCoord();
		this->isAcceptor = true;
	}
	else
	{
		cerr << "not polar atom " << uniqueName << endl;
		exit(0);
	}

	if(atomName.at(0) == 'O')
		this->vdwRadius = 1.6;
	else
		this->vdwRadius = 1.7;

}

PolarAtom::PolarAtom(Residue* res, string atomName, Atom* preC)
{
	if(atomName != "N" || res->triName == "PRO")
	{
		cerr << "constructor for only backbone polar N" << res->triName << "-" <<  atomName << endl;
		exit(0);
	}
	this->uniqueName = res->triName + "-N";
	XYZ n = res->getAtom("N")->getCoord();
	XYZ sup1 = res->getAtom("CA")->getCoord();
	XYZ sup2 = preC->getCoord();
	this->core = n;
	this->support = sup1 + sup2 -n;
	this->isDonor = true;
	this->isAcceptor = false;
	this->vdwRadius = 1.7;
//	this->HCore = true;
}

PolarAtom::~PolarAtom() {
	// TODO Auto-generated destructor stub
}

ProteinChain::ProteinChain() {
	this->pdbID = "xxxx";
	this->chainID = '-';
	this->chainLen = 0;
}

ProteinChain::ProteinChain(char chainID) {
	this->pdbID = "xxxx";
	this->chainID = chainID;
	this->chainLen = 0;
}

ProteinChain::ProteinChain(string pdbID, char chainID){
	this->pdbID = pdbID;
	this->chainID = chainID;
	this->chainLen = 0;
}

void ProteinChain::setPDBID(string pdbID){
	this->pdbID = pdbID;
}

void ProteinChain::setChainID(char c){
	this->chainID = c;
}

string ProteinChain::getPDBID() const{
	return this->pdbID;
}

char ProteinChain::getChainID() const{
	return this->chainID;
}

int ProteinChain::getChainLength() const{
	return this->chainLen;
}

vector<Residue*>& ProteinChain::getResList(){
	return this->resList;
}

Residue* ProteinChain::getResidue(const string& resID) {
	map<string,Residue*>::const_iterator it = resMap.find(resID);
	if(it != resMap.end())
		return it->second;
	else
		return NULL;
}

void ProteinChain::addResidue(Residue* res){
	this->resList.push_back(res);
	this->resMap[res->resID] = res;
	res->setResSeqID(this->resList.size()-1);
	this->chainLen ++;
}

string ProteinChain::getSequence() const{
	string s = "";
	ResName rn = ResName();
	for(int i=0;i<chainLen;i++){
		char c = rn.triToSin(this->resList.at(i)->triName);
		s = s + c;
	}
	return s;
}

int ProteinChain::printPDBFormat(ofstream& out, int startAtomID) const{

	for(int i=0;i<resList.size();i++){
		startAtomID = resList.at(i)->printPDBFormat(out,startAtomID);
	}
	return startAtomID;
}

int ProteinChain::printPDBFormatNoHydrogen(ofstream& out, int startAtomID) const{

	return 0;
}


ProteinChain::~ProteinChain() {
	// TODO Auto-generated destructor stub
}

PDB::PDB() {
	this->pdbID = "pdbx";
	// TODO Auto-generated constructor stub
}

PDB::PDB(const string& pdbFile, const string& pdbID)
{
	this->pdbID = pdbID;
	ifstream input;
	input.open(pdbFile.c_str(),ios::in);

	if (! input.is_open())
    {
        cout << "fail to open file " << pdbFile << endl;
        exit (1);
    }
	init(input);
/*    string s;
    int models = 0;
    int len;
    char lastChainID = '@';
    string lastResID = "XXX";

    char curChainID;
    string curResID;
    char altLoc;
    string resName;

    ProteinChain* curChain;
    Residue* curResidue;

    ResName rn = ResName();
	while(getline(input,s))
    {
        len = s.length();
        if(len < 6) continue;
        string prefix = s.substr(0,6);
        if(prefix == "MODEL ")
        {
            if(models > 0)
                break;
            else
                models = 1;
        }
        if(len < 54) continue;
        if(prefix != "ATOM  " && prefix != "HETATM") continue;
        resName = s.substr(17,3);
        if(resName == "HOH") continue;

        if(!rn.isStandardAminoAcid(resName) && resName != "MSE") continue;


        curChainID = s.at(21);
        curResID = trimString(s.substr(22,5));
        altLoc = s.at(16);
        if(altLoc != ' ' && altLoc != 'A' && altLoc != '1') continue;

        Atom* a = new Atom(s);
        if(curChainID != lastChainID)
        {
             if(getChain(curChainID) == NULL)
             {
                  curChain = new ProteinChain(curChainID);
                  this->chains.push_back(curChain);
             }
             else
            	  curChain = getChain(curChainID);

                  lastChainID = curChainID;
                  lastResID = "XXX";
        }

        if(curResID != lastResID)
        {
        	curResidue = new Residue(curResID,curChainID,resName);
            curChain->addResidue(curResidue);
            residues.push_back(curResidue);
            lastResID = curResID;
        }

        curResidue->addAtom(a);
        curResidue->setAltLoc(altLoc);
    }
*/
    input.close();
}

void PDB::init(std::istream & input)
{
    string s;
    int models = 0;
    int len;
    char lastChainID = '@';
    string lastResID = "XXX";

    char curChainID;
    string curResID;
    char altLoc;
    string resName;

    ProteinChain* curChain;
    Residue* curResidue;

    ResName rn = ResName();
	while(getline(input,s))
    {
        len = s.length();
        if(len < 6) continue;
        string prefix = s.substr(0,6);
        if(prefix == "MODEL ")
        {
            if(models > 0)
                break;
            else
                models = 1;
        }
        if(len < 54) continue;
        if(prefix != "ATOM  " && prefix != "HETATM") continue;
        resName = s.substr(17,3);
        if(resName == "HOH") continue;

        if(!rn.isStandardAminoAcid(resName) && resName != "MSE") continue;


        curChainID = s.at(21);
        curResID = trimString(s.substr(22,5));
        altLoc = s.at(16);
        if(altLoc != ' ' && altLoc != 'A' && altLoc != '1') continue;

        Atom* a = new Atom(s);
        if(curChainID != lastChainID)
        {
             if(getChain(curChainID) == NULL)
             {
                  curChain = new ProteinChain(curChainID);
                  this->chains.push_back(curChain);
             }
             else
            	  curChain = getChain(curChainID);

                  lastChainID = curChainID;
                  lastResID = "XXX";
        }

        if(curResID != lastResID)
        {
        	curResidue = new Residue(curResID,curChainID,resName);
            curChain->addResidue(curResidue);
            residues.push_back(curResidue);
            lastResID = curResID;
        }

        curResidue->addAtom(a);
        curResidue->setAltLoc(altLoc);
    }
}
PDB& PDB::operator=(const PDB& other){
	cout << "operator '=' is inhibited in class PDB" << endl;
	abort();
	return *this;
}

vector<ProteinChain*>& PDB::getChains()
{
    return chains;
}

ProteinChain* PDB::getFirstChain()
{
    if(this->chains.size() > 0)
        return this->chains.at(0);
    return NULL;
}

ProteinChain* PDB::getChain(char c)
{
    for(unsigned int i=0;i<this->chains.size();i++)
    {
        if(this->chains.at(i)->getChainID() == c)
            return chains.at(i);
    }
    return NULL;
}

string PDB::getFirstSeq()
{
    return getFirstChain()->getSequence();
}

vector<Residue*>& PDB::getResList()
{
    return residues;
}

void PDB::printPDBFormat(ofstream& out) const
{
    int startID = 1;
    for(unsigned int i=0;i<this->chains.size();i++)
    {
        ProteinChain* pc = this->chains.at(i);
        for(int j=0;j<pc->getChainLength();j++)
        {
            Residue* res = pc->getResList().at(j);
            startID = res->printPDBFormat(out,startID);
        }
    }
}


PDB::~PDB() {

	// TODO Auto-generated destructor stub
	unsigned int i,j;
	ProteinChain* p;
	for(i = 0;i<chains.size();i++)
	{
		p = chains.at(i);
		delete p;

	}

	Residue* p2;
	Atom* p3;
	for(i=0;i<residues.size();i++)
	{

		p2 = residues.at(i);

		for(j=0;j<p2->getAtomList()->size();j++)
		{
			p3 = p2->getAtomList()->at(j);
			delete p3;

		}
		delete p2;

	}

}

float Phipsi::distance(const Phipsi& other) const{
	float delX = this->phi - other.phi;
	float delY = this->psi - other.psi;
	if(delX > 180)
		delX = 360 - delX;
	else if(delX < -180)
		delX = 360 + delX;

	if(delY > 180)
		delY = 360 - delY;
	else if(delY < -180)
		delY = 360 + delY;
	return sqrt(delX*delX+delY*delY);
}

char Phipsi::regionAB() const{
	if(phi < 0 && psi > -100 && psi < 60)
	{
		return 'A';
	}
	return 'B';
}

Phipsi::~Phipsi() {
	// TODO Auto-generated destructor stub
}

PhipsiLib::PhipsiLib(){
	this->pointNum = 200;
	this->ppList.push_back(new Phipsi( -62.71 ,  -41.57));
	this->ppList.push_back(new Phipsi( -60.44 ,  -36.89));
	this->ppList.push_back(new Phipsi( -65.79 ,  -36.65));
	this->ppList.push_back(new Phipsi( -67.17 ,  -41.65));
	this->ppList.push_back(new Phipsi( -59.91 ,  -45.30));
	this->ppList.push_back(new Phipsi( -56.70 ,  -41.11));
	this->ppList.push_back(new Phipsi( -64.61 ,  -46.65));
	this->ppList.push_back(new Phipsi( -59.68 ,  -50.44));
	this->ppList.push_back(new Phipsi( -54.61 ,  -46.94));
	this->ppList.push_back(new Phipsi( -71.59 ,  -46.40));
	this->ppList.push_back(new Phipsi( -72.03 ,  -37.93));
	this->ppList.push_back(new Phipsi( -62.47 ,  -30.49));
	this->ppList.push_back(new Phipsi( -54.86 ,  -32.60));
	this->ppList.push_back(new Phipsi( -69.62 ,  -31.26));
	this->ppList.push_back(new Phipsi( -48.60 ,  -39.19));
	this->ppList.push_back(new Phipsi( -65.64 ,  -55.36));
	this->ppList.push_back(new Phipsi( -53.82 ,  -55.75));
	this->ppList.push_back(new Phipsi( -46.67 ,  -48.77));
	this->ppList.push_back(new Phipsi( -57.50 ,  -24.42));
	this->ppList.push_back(new Phipsi( -66.32 ,  -23.90));
	this->ppList.push_back(new Phipsi( -78.56 ,  -32.35));
	this->ppList.push_back(new Phipsi( -81.49 ,  -42.96));
	this->ppList.push_back(new Phipsi( -75.23 ,  -23.39));
	this->ppList.push_back(new Phipsi( -62.73 ,  -17.29));
	this->ppList.push_back(new Phipsi( -70.93 ,  -15.55));
	this->ppList.push_back(new Phipsi( -80.44 ,  -56.71));
	this->ppList.push_back(new Phipsi( -86.95 ,  -23.82));
	this->ppList.push_back(new Phipsi( -92.42 ,  -35.65));
	this->ppList.push_back(new Phipsi( -81.19 ,  -14.33));
	this->ppList.push_back(new Phipsi( -66.67 ,   -7.74));
	this->ppList.push_back(new Phipsi( -76.53 ,   -6.37));
	this->ppList.push_back(new Phipsi( -96.72 ,  -50.52));
	this->ppList.push_back(new Phipsi( -91.95 ,  -12.06));
	this->ppList.push_back(new Phipsi(-100.84 ,  -23.37));
	this->ppList.push_back(new Phipsi( -85.93 ,   -3.37));
	this->ppList.push_back(new Phipsi(-108.71 ,  -38.01));
	this->ppList.push_back(new Phipsi( -76.47 ,    6.23));
	this->ppList.push_back(new Phipsi(-104.47 ,  -10.15));
	this->ppList.push_back(new Phipsi( -96.36 ,   -0.46));
	this->ppList.push_back(new Phipsi( -89.35 ,    6.81));
	this->ppList.push_back(new Phipsi(-116.16 ,  -21.63));
	this->ppList.push_back(new Phipsi( -38.10 ,  -59.62));
	this->ppList.push_back(new Phipsi(-107.80 ,    3.09));
	this->ppList.push_back(new Phipsi(-100.12 ,   11.45));
	this->ppList.push_back(new Phipsi(-117.52 ,   -6.73));
	this->ppList.push_back(new Phipsi(-117.12 ,  -57.58));
	this->ppList.push_back(new Phipsi( -95.66 ,  -77.23));
	this->ppList.push_back(new Phipsi( -91.15 ,   20.75));
	this->ppList.push_back(new Phipsi(-129.66 ,  -36.11));
	this->ppList.push_back(new Phipsi(-112.26 ,   14.81));
	this->ppList.push_back(new Phipsi(-133.22 ,  -11.58));
	this->ppList.push_back(new Phipsi(-124.53 ,    6.65));
	this->ppList.push_back(new Phipsi(-105.87 ,   25.91));
	this->ppList.push_back(new Phipsi( -58.02 ,  -83.01));
	this->ppList.push_back(new Phipsi(-124.36 ,   23.11));
	this->ppList.push_back(new Phipsi( -78.31 ,   41.41));
	this->ppList.push_back(new Phipsi(-142.95 ,  -64.95));
	this->ppList.push_back(new Phipsi(-143.81 ,   13.38));
	this->ppList.push_back(new Phipsi(-119.01 ,   38.30));
	this->ppList.push_back(new Phipsi(-119.93 ,  -94.85));
	this->ppList.push_back(new Phipsi(-138.92 ,   38.56));
	this->ppList.push_back(new Phipsi(-172.29 ,  -27.81));
	this->ppList.push_back(new Phipsi(  -5.37 ,  -81.80));
	this->ppList.push_back(new Phipsi(-109.17 ,  132.79));
	this->ppList.push_back(new Phipsi(-108.81 ,  123.49));
	this->ppList.push_back(new Phipsi(-117.57 ,  127.44));
	this->ppList.push_back(new Phipsi(-118.05 ,  136.97));
	this->ppList.push_back(new Phipsi(-118.00 ,  118.36));
	this->ppList.push_back(new Phipsi(-125.28 ,  131.48));
	this->ppList.push_back(new Phipsi(-127.56 ,  122.79));
	this->ppList.push_back(new Phipsi(-108.62 ,  114.24));
	this->ppList.push_back(new Phipsi(-100.17 ,  128.43));
	this->ppList.push_back(new Phipsi( -99.45 ,  118.02));
	this->ppList.push_back(new Phipsi(-110.11 ,  143.24));
	this->ppList.push_back(new Phipsi(-100.40 ,  139.45));
	this->ppList.push_back(new Phipsi(-128.65 ,  139.96));
	this->ppList.push_back(new Phipsi(-121.40 ,  146.83));
	this->ppList.push_back(new Phipsi(-134.95 ,  131.25));
	this->ppList.push_back(new Phipsi(-127.40 ,  111.21));
	this->ppList.push_back(new Phipsi(-115.32 ,  105.65));
	this->ppList.push_back(new Phipsi(-138.60 ,  119.63));
	this->ppList.push_back(new Phipsi(-138.74 ,  141.26));
	this->ppList.push_back(new Phipsi(-131.24 ,  150.02));
	this->ppList.push_back(new Phipsi(-113.22 ,  153.33));
	this->ppList.push_back(new Phipsi( -91.01 ,  134.29));
	this->ppList.push_back(new Phipsi( -90.81 ,  123.32));
	this->ppList.push_back(new Phipsi(-101.02 ,  152.94));
	this->ppList.push_back(new Phipsi(-123.72 ,  157.97));
	this->ppList.push_back(new Phipsi( -90.39 ,  145.82));
	this->ppList.push_back(new Phipsi(-100.97 ,  105.75));
	this->ppList.push_back(new Phipsi( -90.37 ,  110.74));
	this->ppList.push_back(new Phipsi( -81.99 ,  129.04));
	this->ppList.push_back(new Phipsi( -80.58 ,  139.50));
	this->ppList.push_back(new Phipsi( -81.23 ,  117.30));
	this->ppList.push_back(new Phipsi(-110.58 ,  164.44));
	this->ppList.push_back(new Phipsi( -88.14 ,  157.53));
	this->ppList.push_back(new Phipsi( -78.65 ,  149.95));
	this->ppList.push_back(new Phipsi( -73.22 ,  125.73));
	this->ppList.push_back(new Phipsi( -96.67 ,  167.33));
	this->ppList.push_back(new Phipsi( -71.22 ,  135.23));
	this->ppList.push_back(new Phipsi( -70.37 ,  144.53));
	this->ppList.push_back(new Phipsi(-134.04 ,  160.70));
	this->ppList.push_back(new Phipsi(-140.74 ,  152.22));
	this->ppList.push_back(new Phipsi(-123.92 ,  169.98));
	this->ppList.push_back(new Phipsi( -76.83 ,  160.13));
	this->ppList.push_back(new Phipsi(-146.93 ,  130.97));
	this->ppList.push_back(new Phipsi( -81.98 ,  169.43));
	this->ppList.push_back(new Phipsi( -68.82 ,  154.34));
	this->ppList.push_back(new Phipsi(-108.56 , -179.84));
	this->ppList.push_back(new Phipsi(-149.54 ,  142.95));
	this->ppList.push_back(new Phipsi(-136.72 ,  172.17));
	this->ppList.push_back(new Phipsi(-145.19 ,  163.00));
	this->ppList.push_back(new Phipsi(-151.06 ,  154.55));
	this->ppList.push_back(new Phipsi(-108.83 ,   93.20));
	this->ppList.push_back(new Phipsi( -68.43 ,  116.31));
	this->ppList.push_back(new Phipsi( -78.79 ,  103.37));
	this->ppList.push_back(new Phipsi( -63.16 ,  128.77));
	this->ppList.push_back(new Phipsi( -90.94 ,   95.27));
	this->ppList.push_back(new Phipsi( -61.74 ,  138.40));
	this->ppList.push_back(new Phipsi( -61.79 ,  147.62));
	this->ppList.push_back(new Phipsi( -68.50 ,  165.94));
	this->ppList.push_back(new Phipsi(-126.07 ,   94.90));
	this->ppList.push_back(new Phipsi( -90.20 , -176.79));
	this->ppList.push_back(new Phipsi(-141.05 ,  103.07));
	this->ppList.push_back(new Phipsi( -58.84 ,  157.36));
	this->ppList.push_back(new Phipsi( -53.54 ,  132.85));
	this->ppList.push_back(new Phipsi( -52.94 ,  143.20));
	this->ppList.push_back(new Phipsi( -74.20 ,  179.78));
	this->ppList.push_back(new Phipsi( -54.33 ,  121.37));
	this->ppList.push_back(new Phipsi(-126.43 , -170.97));
	this->ppList.push_back(new Phipsi(-155.86 ,  116.09));
	this->ppList.push_back(new Phipsi(-160.77 ,  148.35));
	this->ppList.push_back(new Phipsi(-162.21 ,  134.10));
	this->ppList.push_back(new Phipsi(-153.57 ,  171.62));
	this->ppList.push_back(new Phipsi(-159.20 ,  161.49));
	this->ppList.push_back(new Phipsi(-145.25 , -175.88));
	this->ppList.push_back(new Phipsi( -80.19 ,   82.18));
	this->ppList.push_back(new Phipsi(-122.33 ,   77.75));
	this->ppList.push_back(new Phipsi( -93.64 ,   75.65));
	this->ppList.push_back(new Phipsi(-138.24 ,   79.85));
	this->ppList.push_back(new Phipsi( -54.25 ,   96.28));
	this->ppList.push_back(new Phipsi( -40.61 ,  127.84));
	this->ppList.push_back(new Phipsi(-104.40 , -156.29));
	this->ppList.push_back(new Phipsi( -81.68 , -159.28));
	this->ppList.push_back(new Phipsi(-173.28 ,  156.71));
	this->ppList.push_back(new Phipsi(-167.83 ,  170.51));
	this->ppList.push_back(new Phipsi(-159.27 ,   91.56));
	this->ppList.push_back(new Phipsi(-162.77 , -174.54));
	this->ppList.push_back(new Phipsi(-148.38 , -155.90));
	this->ppList.push_back(new Phipsi( -81.65 ,   64.08));
	this->ppList.push_back(new Phipsi(-130.51 ,   61.59));
	this->ppList.push_back(new Phipsi(-152.79 ,   64.32));
	this->ppList.push_back(new Phipsi(-100.02 ,   51.12));
	this->ppList.push_back(new Phipsi(-122.55 , -135.41));
	this->ppList.push_back(new Phipsi( -92.38 , -124.40));
	this->ppList.push_back(new Phipsi(-157.80 , -114.40));
	this->ppList.push_back(new Phipsi(  79.97 ,   11.05));
	this->ppList.push_back(new Phipsi(  67.21 ,   13.11));
	this->ppList.push_back(new Phipsi(  75.20 ,   24.73));
	this->ppList.push_back(new Phipsi(  77.13 ,   -1.28));
	this->ppList.push_back(new Phipsi(  92.68 ,   16.58));
	this->ppList.push_back(new Phipsi(  91.44 ,    1.17));
	this->ppList.push_back(new Phipsi(  60.91 ,   25.89));
	this->ppList.push_back(new Phipsi(  67.07 ,   37.97));
	this->ppList.push_back(new Phipsi(  55.06 ,   37.75));
	this->ppList.push_back(new Phipsi(  87.84 ,  -12.94));
	this->ppList.push_back(new Phipsi( 111.79 ,    6.45));
	this->ppList.push_back(new Phipsi( 102.18 ,  -12.68));
	this->ppList.push_back(new Phipsi(  97.17 ,   46.11));
	this->ppList.push_back(new Phipsi(  58.70 ,   52.03));
	this->ppList.push_back(new Phipsi(  48.74 ,   46.96));
	this->ppList.push_back(new Phipsi(  43.41 ,   60.83));
	this->ppList.push_back(new Phipsi(  62.67 ,   78.25));
	this->ppList.push_back(new Phipsi( 104.48 ,  -30.61));
	this->ppList.push_back(new Phipsi(  73.14 ,  -44.39));
	this->ppList.push_back(new Phipsi( 131.11 ,  -17.09));
	this->ppList.push_back(new Phipsi( 160.05 ,   42.52));
	this->ppList.push_back(new Phipsi( 142.52 ,  -65.36));
	this->ppList.push_back(new Phipsi(   4.13 ,   98.81));
	this->ppList.push_back(new Phipsi( 123.61 , -160.45));
	this->ppList.push_back(new Phipsi(  92.55 , -154.26));
	this->ppList.push_back(new Phipsi( 100.14 ,  178.37));
	this->ppList.push_back(new Phipsi( 128.04 ,  168.58));
	this->ppList.push_back(new Phipsi(  76.51 , -174.33));
	this->ppList.push_back(new Phipsi(  84.21 ,  162.79));
	this->ppList.push_back(new Phipsi( 103.76 ,  146.02));
	this->ppList.push_back(new Phipsi(  67.39 , -154.11));
	this->ppList.push_back(new Phipsi( 152.15 , -168.57));
	this->ppList.push_back(new Phipsi( 114.53 , -122.60));
	this->ppList.push_back(new Phipsi( 146.81 , -140.10));
	this->ppList.push_back(new Phipsi(  77.05 , -123.32));
	this->ppList.push_back(new Phipsi(  56.48 , -134.67));
	this->ppList.push_back(new Phipsi( 159.19 ,  165.30));
	this->ppList.push_back(new Phipsi( 175.66 , -154.48));
	this->ppList.push_back(new Phipsi( 176.13 , -177.66));
	this->ppList.push_back(new Phipsi( 147.35 ,  128.99));
	this->ppList.push_back(new Phipsi(  52.69 , -117.81));
	this->ppList.push_back(new Phipsi( 106.36 ,  108.41));
	this->ppList.push_back(new Phipsi(  68.96 ,  118.77));
	this->ppList.push_back(new Phipsi(  81.03 ,  -73.88));
	creatIndexTable();
}
int PhipsiLib::findNearestPointsWithoutIndex(Phipsi* pp)
{
	float minDis = 10000.0;
	int minIndex = -1;
	float d;
	for(int i=0;i<pointNum;i++)
	{
		d = pp->distance(*(ppList.at(i)));
		//cout << d << endl;
		if(d < minDis)
		{
			minDis = d;
			minIndex = i;
		}
	}
	return minIndex;
}

void PhipsiLib::creatIndexTable()
{
	for(int i=0;i<36;i++)
	{
		for(int j=0;j<36;j++)
		{
			for(int k=0;k<20;k++)
			{
				this->indexTable[i][j][k] = -1;
			}
		}
	}

	set<unsigned int> possibleNeighbor;
	Phipsi *pp, *pp1, *pp2, *pp3, *pp4;
	set<unsigned int>::const_iterator it;

	for(int i=0;i<36;i++)
	{
		for(int j=0;j<36;j++)
		{
			possibleNeighbor.clear();
			float phi0 = i*10-180;
			float psi0 = j*10-180;
			for(unsigned int a=0;a<ppList.size();a++)
			{
				pp = ppList.at(a);
				if((pp->phi) > phi0 && (pp->phi) < phi0+10 && (pp->psi) > psi0 && (pp->psi) < psi0+10)
					possibleNeighbor.insert(a);
			}

			for(int a=0;a<11;a++)
			{
				pp1 = new Phipsi(phi0+a,psi0);
				pp2 = new Phipsi(phi0+a, psi0+10);
				pp3 = new Phipsi(phi0, psi0+a);
				pp4 = new Phipsi(phi0+10, psi0+a);
				possibleNeighbor.insert(findNearestPointsWithoutIndex(pp1));
				possibleNeighbor.insert(findNearestPointsWithoutIndex(pp2));
				possibleNeighbor.insert(findNearestPointsWithoutIndex(pp3));
				possibleNeighbor.insert(findNearestPointsWithoutIndex(pp4));
				delete pp1;
				delete pp2;
				delete pp3;
				delete pp4;
			}

			int n=0;
			for(it = possibleNeighbor.begin();it!=possibleNeighbor.end();it++)
			{
				int k = *it;
				this->indexTable[i][j][n] = k;
				n++;
			}
		}
	}
}

Phipsi* PhipsiLib::indexToPhipsi(int id) const
{
	return this->ppList.at(id);
}

int PhipsiLib::phipsiToIndex(const Phipsi* pp) const
{
	int i = (int)((pp->phi+180)/10);
	int j = (int)((pp->psi+180)/10);
	double minD = 10000.0;
	int minIndex = -1;
	for(int k=0;k<20;k++)
	{
		int id = this->indexTable[i][j][k];
		if(id < 0)
			break;
		double dist = pp->distance(* ppList.at(id));
		if(dist < minD)
		{
			minD = dist;
			minIndex = id;
		}
	}
	return minIndex;
}

PhipsiLib::~PhipsiLib() {
	// TODO Auto-generated destructor stub
	for(int i=0;i<this->pointNum;i++)
	{
		delete this->ppList.at(i);
	}
}

ResPairOrientation::ResPairOrientation(){
	for(int i=0;i<10;i++){
		points.push_back(XYZ(0,0,0));
	}
}

ResPairOrientation::ResPairOrientation(LocalFrame& csA, LocalFrame& csB){

	LocalFrame newCsB = buildLocalFrame(csA, csB);
	addAtoms(newCsB);
}

ResPairOrientation::ResPairOrientation(LocalFrame& csB){
	addAtoms(csB);
}

ResPairOrientation::ResPairOrientation(NSPproteinrep::BackBoneSite& resA, NSPproteinrep::BackBoneSite& resB){
	LocalFrame csA = getBackboneSiteLocalFrame(resA);
	LocalFrame csB = getBackboneSiteLocalFrame(resB);
	LocalFrame newCsB = buildLocalFrame(csA, csB);
	addAtoms(newCsB);
}

ResPairOrientation::ResPairOrientation(string& s){
	vector<string> spt;
	splitString(s," ",&spt);
	if(spt.size() != 15){
		cout << "orientation string error: " + s << endl;
		exit(1);
	}

	points.push_back(XYZ(0.826, -1.204, 0));
	points.push_back(XYZ(0,0,0));
	points.push_back(XYZ(0.862, 1.257, 0));
	points.push_back(XYZ(-1.040, 0.006, 1.345));
	points.push_back(XYZ(-2.332, -0.663, 0.981));
	points.push_back(XYZ(atof(spt[0].c_str()), atof(spt[1].c_str()), atof(spt[2].c_str())));
	points.push_back(XYZ(atof(spt[3].c_str()), atof(spt[4].c_str()), atof(spt[5].c_str())));
	points.push_back(XYZ(atof(spt[6].c_str()), atof(spt[7].c_str()), atof(spt[8].c_str())));
	points.push_back(XYZ(atof(spt[9].c_str()), atof(spt[10].c_str()), atof(spt[11].c_str())));
	points.push_back(XYZ(atof(spt[12].c_str()), atof(spt[13].c_str()), atof(spt[14].c_str())));

}


void ResPairOrientation::addAtoms(LocalFrame& csB){

	points.push_back(XYZ(0.826, -1.204, 0));
	points.push_back(XYZ(0,0,0));
	points.push_back(XYZ(0.862, 1.257, 0));
	points.push_back(XYZ(-1.040, 0.006, 1.345));
	points.push_back(XYZ(-2.332, -0.663, 0.981));

	points.push_back(csB.local2globalcrd(points[0]));
	points.push_back(csB.local2globalcrd(points[1]));
	points.push_back(csB.local2globalcrd(points[2]));
	points.push_back(csB.local2globalcrd(points[3]));
	points.push_back(csB.local2globalcrd(points[4]));

}

double ResPairOrientation::rmsd(const ResPairOrientation& other) const{
	QuatFit qf;
	vector<double> weights;
	for(int i=0;i<10;i++)
		weights.push_back(1.0);

	return qf.setup(this->points, other.points, weights);
}

XYZ ResPairOrientation::getTern(int index){
	if(index >=0 && index < 10)
		return points[index];
	else {
		cout << "index out of range" << endl;
		exit(1);
	}
}

string ResPairOrientation::toString() const{
	string rpoString = "";
	char s[30];
	for(int i=0;i<10;i++){
		sprintf(s, " %7.3f %7.3f %7.3f", points[i][0], points[i][1], points[i][2]);
		string ss(s);
		rpoString = rpoString + ss;
	}
	return rpoString;
}

BackboneSitesPair::BackboneSitesPair(BackBoneSite* siteA, BackBoneSite* siteB){
	if(siteA->resseq < siteB->resseq)
	{
		this->siteA = siteA;
		this->siteB = siteB;
	}
	else
	{
		this->siteA = siteB;
		this->siteB = siteA;
	}
	int sep = this->siteB->resid - this->siteA->resid;
	if(sep > 5) sep = 5;
	this->seqSep = sep;

	this->ori = ResPairOrientation(*siteA, *siteB);
}




LocalFrame getBackboneSiteLocalFrame(const NSPproteinrep::BackBoneSite& bbSite){
	XYZ n = bbSite.ncrd();
	XYZ ca = bbSite.cacrd();
	XYZ c = bbSite.ccrd();

	XYZ can = ~(n-ca);
	XYZ cac = ~(c-ca);
	XYZ z = ~(can^cac);
	XYZ x = ~(can+cac);
	XYZ y = ~(z^x);

	LocalFrame cs;
	cs.origin_ = ca;
	cs.axis_.push_back(x);
	cs.axis_.push_back(y);
	cs.axis_.push_back(z);
	return cs;
}

LocalFrame buildLocalFrame(LocalFrame& csNewGlobal, LocalFrame& csOld){
	LocalFrame csNew;
    csNew.origin_ = csNewGlobal.global2localcrd(csOld.origin_);
    double c[3][3];
    double bColumnJ[3];
    double aRowI[3];
    for(int j=0;j<3;j++){
    	for(int k=0;k<3;k++){
    		bColumnJ[k] = csOld.axis_[j][k];
    	}
    	for(int i=0;i<3;i++){
    		aRowI[0] = csNewGlobal.axis_[i][0];
    		aRowI[1] = csNewGlobal.axis_[i][1];
    		aRowI[2] = csNewGlobal.axis_[i][2];
    		double s = 0;
    		for(int k=0;k<3;k++){
    			s += aRowI[k]*bColumnJ[k];
    		}
    		c[i][j] = s;
    	}
    }
    csNew.axis_.push_back(XYZ(c[0][0], c[1][0], c[2][0]));
    csNew.axis_.push_back(XYZ(c[0][1], c[1][1], c[2][1]));
    csNew.axis_.push_back(XYZ(c[0][2], c[1][2], c[2][2]));
    return csNew;
}


} /* namespace NSPdesignseq */
