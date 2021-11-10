/*
 * StructureInfo.cpp
 *
 *  Created on: 2017��11��1��
 *      Author: notxp
 */

#include "designseq/StructureInfo.h"

namespace NSPdesignseq {

BackboneHBond::BackboneHBond(XYZ& N, XYZ& H, XYZ& O, Residue* resA, Residue* resB) {
	this->donerChainID = resA->chainID;
	this->acceptorChainID = resB->chainID;
	this->donerResID = resA->resID;
	this->acceptorResID = resB->resID;
	if(donerChainID != acceptorChainID)
		this->seqSeparation = 999;
	else
		this->seqSeparation = abs(resA->resSeqID - resB->resSeqID);

	this->hDoner = N;
	this->hydrogen = H;
	this->hAcceptor = O;
	this->distance = hydrogen.distance(hAcceptor);
	this->angle = NSPdesignseq::angle(N,H,O);
}

string BackboneHBond::toString()
{
	char s[100];
	sprintf(s,"%-4d %c%-4s %c%-4s %5.3f %6.2f", this->seqSeparation, donerChainID, donerResID.c_str(), acceptorChainID, acceptorResID.c_str(), distance, angle);
	return s;
}

StructureInfo::StructureInfo(PDB* protein) {
//	cout << "start init StructureInfo: " << endl;
	vector<Residue*> tmpResList = protein->getResList();
	int index = 0;
	char s[10];
	//cout << "tmpList " << tmpResList.size() << endl;
	for(unsigned int i=0;i<tmpResList.size();i++)
	{
		Residue* res = tmpResList.at(i);
		if(res->hasThreeCoreAtoms() && res->hasAtom("O"))
		{
			resList.push_back(res);
			NList.push_back(&res->getAtom("N")->getCoord());
			CAList.push_back(&res->getAtom("CA")->getCoord());
			CList.push_back(&res->getAtom("C")->getCoord());
			sprintf(s,"%c%s",res->chainID,res->resID.c_str());
			chainIDResIDToSeqID[string(s)] = index;
			index++;
		}
	}
//	cout << "res num: " << resList.size() << endl;

	this->resNum = resList.size();
	this->ssSeq = new char[resNum+1];
	for(int i=0;i<resNum;i++)
	{
		ssSeq[i] = 'C';
	}
	ssSeq[resNum] = '\0';

	phiList = new float[resNum];
	psiList = new float[resNum];
	omgList = new float[resNum];

//	cout << "finish init" << endl;
}

StructureInfo::StructureInfo(ProteinChain* chain) {
	vector<Residue*> tmpResList = chain->getResList();
	int index = 0;
	char s[10];
	//cout << "tmpList " << tmpResList.size() << endl;
	for(unsigned int i=0;i<tmpResList.size();i++)
	{
		Residue* res = tmpResList.at(i);
		if(res->hasThreeCoreAtoms()&& res->hasAtom("O"))
		{
			resList.push_back(res);
			NList.push_back(&res->getAtom("N")->getCoord());
			CAList.push_back(&res->getAtom("CA")->getCoord());
			CList.push_back(&res->getAtom("C")->getCoord());
			sprintf(s,"%c%s",res->chainID,res->resID.c_str());
			chainIDResIDToSeqID[string(s)] = index;
			index++;
		}
	}
	this->resNum = resList.size();
	this->ssSeq = new char[resNum+1];
	for(int i=0;i<resNum;i++)
	{
		ssSeq[i] = 'C';
	}
	ssSeq[resNum] = '\0';

	phiList = new float[resNum];
	psiList = new float[resNum];
	omgList = new float[resNum];

}

StructureInfo::StructureInfo(vector<Residue*>& residueList) {

	int index = 0;
	char s[10];
	//cout << "tmpList " << tmpResList.size() << endl;
	for(unsigned int i=0;i<residueList.size();i++)
	{
		Residue* res = residueList.at(i);
		if(res->hasThreeCoreAtoms()&& res->hasAtom("O"))
		{
			resList.push_back(res);
			NList.push_back(&res->getAtom("N")->getCoord());
			CAList.push_back(&res->getAtom("CA")->getCoord());
			CList.push_back(&res->getAtom("C")->getCoord());
			sprintf(s,"%c%s",res->chainID,res->resID.c_str());
			chainIDResIDToSeqID[string(s)] = index;
			index++;
		}
	}
	this->resNum = resList.size();
	this->ssSeq = new char[resNum+1];
	for(int i=0;i<resNum;i++)
	{
		ssSeq[i] = 'C';
	}
	ssSeq[resNum] = '\0';

	phiList = new float[resNum];
	psiList = new float[resNum];
	omgList = new float[resNum];

}

void StructureInfo::updateTorsion()
{
//	cout << "update torsion:" << endl;
	phiList[0] = 360.0;
	omgList[0] = 360.0;

	/*
	 * first residue
	 */
	if(resNum <2)
	{
		phiList[0] = 360.0;
		return;
	}
	else
	{
		XYZ A = *NList.at(0);
		XYZ B = *CAList.at(0);
		XYZ C = *CList.at(0);
		XYZ D = *NList.at(1);
		if(C.distance(D) < 1.6)
			psiList[0] = dihedral(A,B,C,D);
		else
			psiList[0] = 360.0;
	}

	/*
	 * second residue last but one residue
	 */
	double phi,psi,omg;
	XYZ A,B,C,D,E,F;
	int i;

	for(i=1;i<resNum-1;i++)
	{
		phi=360.0;
		psi=360.0;
		omg=360.0;
		A = *CAList.at(i-1);
		B = *CList.at(i-1);
		C = *NList.at(i);
		D = *CAList.at(i);
		E = *CList.at(i);
		F = *NList.at(i+1);

		if(B.distance(C) < 1.6)
		{
			phi = dihedral(B,C,D,E);
			omg = dihedral(A,B,C,D);
		}
		if(E.distance(F) < 1.6)
		{
			psi = dihedral(C,D,E,F);
		}
		phiList[i] = phi;
		psiList[i] = psi;
		omgList[i] = omg;
	}

	phi = 360.0;
	psi = 360.0;
	omg = 360.0;
	A = *CAList.at(i-1);
	B = *CList.at(i-1);
	C = *NList.at(i);
	D = *CAList.at(i);
	E = *CList.at(i);
	if(B.distance(C) < 1.6)
	{
		phi = dihedral(B,C,D,E);
		omg = dihedral(A,B,C,D);
	}
	phiList[i] = phi;
	psiList[i] = psi;
	omgList[i] = omg;


}

void StructureInfo::updateBBHBonds()
{
	vector<XYZ*> HNList;
	vector<XYZ*> OList;

	for(int i=0;i<resNum;i++)
	{
		Residue* res = resList.at(i);
		Atom* O = res->getAtom("O");
		if(O != NULL)
			OList.push_back(&(O->getCoord()));
		else
			OList.push_back(NULL);

		if(i==0 || res->triName == "PRO")
			HNList.push_back(NULL);
		else
		{
			Atom* aC1 = resList.at(i-1)->getAtom("C");
			Atom* aC2 = resList.at(i)->getAtom("CA");
			if(aC1->getCoord().distance(aC2->getCoord()) > 3.0)
			{
				HNList.push_back(NULL);
			}
			else
			{
				XYZ NC1 = ~(aC1->getCoord() - *NList[i]);
				XYZ NC2 = ~(aC2->getCoord() - *NList[i]);
				XYZ NH = ~(NC1 + NC2);
				XYZ* H = new XYZ(*NList[i] - NH); //need delete
				HNList.push_back(H);

			}
		}
	}

	int i,j;
	XYZ *NI, *NJ, *OI, *OJ, *HI, *HJ;
	float d1,d2;
	for(i=1;i<resNum;i++)
	{
		NI = NList.at(i);
		OI = OList.at(i);
		HI = HNList.at(i);

		for(j=i+2;j<resNum;j++)
		{
			NJ = NList.at(j);
			OJ = OList.at(j);
			HJ = HNList.at(j);




			if(HI != NULL && OJ != NULL)
			{
				d1 = HI->distance(*OJ);
				if(d1 < 2.5 && d1 > 1.6)
				{
					BackboneHBond* hb = new BackboneHBond(*NI,*HI, *OJ, resList.at(i), resList.at(j));
					if(hb->angle > 100)
						hbList.push_back(hb);
					else
						delete hb;
				}
			//	cout << "HI OJ" << HI->distance(*OJ) << endl;

			}
			if(HJ != NULL && OI != NULL)
			{
				d2 = HJ->distance(*OI);
				if(d2 < 2.5 && d2 > 1.6)
				{
					BackboneHBond* hb = new BackboneHBond(*NJ, *HJ, *OI, resList.at(j), resList.at(i));
					if(hb->angle > 100)
						hbList.push_back(hb);
					else
						delete hb;
				}
			//	cout << "HJ OI" << HJ->distance(*OI) << endl;
			}
		}
	}

	XYZ* t;
	for(i=0;i<resNum;i++)
	{
		t = HNList.at(i);
		if(t != NULL)
			delete t;
	}


}

void StructureInfo::updateSecondaryStructure()
{
//	cout << "update ss" << endl;
	updateBBHBonds();

	int NHbList[resNum];
	int OHbList[resNum];
	for(int i=0;i<resNum;i++)
	{
		NHbList[i] = 0;
		OHbList[i] = 0;
	}


	int resASeqID, resBSeqID;
	char s1[10], s2[10];

	for(unsigned int i=0;i<hbList.size();i++)
	{
		BackboneHBond *hb = hbList.at(i);
		sprintf(s1,"%c%s",hb->donerChainID, hb->donerResID.c_str());
		sprintf(s2, "%c%s", hb->acceptorChainID, hb->acceptorResID.c_str());
		resASeqID = chainIDResIDToSeqID[s1];
		resBSeqID = chainIDResIDToSeqID[s2];
		if(NHbList[resASeqID] != 0)
		{
			int oldResBSeqID = resASeqID + NHbList[resASeqID];
			if(abs(resASeqID - resBSeqID) > abs(resASeqID - oldResBSeqID))
				continue;
		}

		if(OHbList[resBSeqID] != 0)
		{
			int oldSeqASeqID =	resBSeqID + OHbList[resBSeqID];
			if(abs(resASeqID - resBSeqID) > abs(oldSeqASeqID - resBSeqID))
				continue;
		}

		NHbList[resASeqID] = resBSeqID - resASeqID;
		OHbList[resBSeqID] = resASeqID - resBSeqID;
	}
/*
	for(int i=0;i<resNum;i++)
	{
		printf("%-3d %4d %4d\n",i+1, NHbList[i], OHbList[i]);
	}
*/

	/*
	 * assign Helix
	 */
	for(int i=1;i<resNum-1;i++)
	{
		if(OHbList[i+1] == 4 && OHbList[i] == 4)
			this->ssSeq[i] = 'H';
		else if(NHbList[i] == -4 && NHbList[i+1] != -4)
			this->ssSeq[i] = 'C';
		else if(this->ssSeq[i-1] == 'H')
			this->ssSeq[i] = 'H';
		else if(OHbList[i-1] == 3 && OHbList[i] == 3)
			this->ssSeq[i] = 'G';
		else if(NHbList[i] == -3 && NHbList[i+1] != -3)
			this->ssSeq[i] = 'C';
		else if(this->ssSeq[i-1] == 'G')
			this->ssSeq[i] = 'G';
	}

	/*
	 * assign beta
	 */
	for(int i=1;i<resNum-1;i++)
	{
		if(NHbList[i] == OHbList[i] && abs(NHbList[i]) > 3)
			this->ssSeq[i] = 'E';
		if(NHbList[i]+2 == OHbList[i] && abs(NHbList[i]) > 3)
			this->ssSeq[i] = 'E';
	}

	for(int i=1;i<resNum-1;i++)
	{
		if(this->ssSeq[i] == 'C' && (this->ssSeq[i-1] == 'E' || this->ssSeq[i+1] == 'E'))
			this->ssSeq[i] = 'e';
	}

	for(int i=1;i<resNum-1;i++)
	{
		if(this->ssSeq[i] == 'e')
			this->ssSeq[i] = 'E';
	}
}

void StructureInfo::updateSAI(SasaPSD* rsp)
{
//	cout << "update sasa" << endl;
	for(int i=0;i<resNum;i++)
	{
		int n = rsp->exposeNum(resList.at(i), resList);
		float sai = rsp->pointNumToIndex(n);
	//	cout << "UPDATESAI: " <<  resList.at(i)->resID << " " << n << " " << sai << endl;

		saiList.push_back(rsp->pointNumToIndex(n));
	}
//	cout << "finish update sasa" << endl;
}

string StructureInfo::getSecondaryStructureSeq()
{
	return string(this->ssSeq);
}

Phipsi StructureInfo::getPhipsi(Residue* res)
{
	char s[10];
	sprintf(s,"%c%s",res->chainID,res->resID.c_str());
	map<string,int>::const_iterator it;
	it = chainIDResIDToSeqID.find(string(s));
	if(it != chainIDResIDToSeqID.end())
	{
		int id = it->second;
		Phipsi pp(phiList[id],psiList[id]);
		return Phipsi(phiList[id],psiList[id]);
	}
	else
	{
		cerr << "can't find phi-psi of Residue: " << s << endl;
		return Phipsi();
	}
}

float StructureInfo::getSai(Residue* res)
{
	char s[10];
	sprintf(s,"%c%s",res->chainID,res->resID.c_str());
	map<string,int>::const_iterator it;
	it = chainIDResIDToSeqID.find(string(s));
	if(it != chainIDResIDToSeqID.end())
	{
		int id = it->second;
		return this->saiList.at(id);
	}
	else
	{
		cerr << "can't find sai of Residue: " << s << endl;
		return 0.5;
	}
}

StructureInfo::~StructureInfo() {
	// TODO Auto-generated destructor stub

	for(unsigned int i=0;i<hbList.size();i++)
	{
		delete hbList.at(i);
	}

	delete [] ssSeq;
	delete [] phiList;
	delete [] psiList;
	delete [] omgList;
}

void backboneSiteListUpdateSasa(vector<BackBoneSite*>& bsList){
	vector<Residue*> resList;
	char s[20];
	for(int i=0;i<bsList.size();i++){
		BackBoneSite* bs = bsList.at(i);
		sprintf(s,"%d",i+1);
		string resid = string(s);
		Residue* res = new Residue(resid, 'A', bs->resname);
		res->addAtom(new Atom("N", bs->ncrd()));
		res->addAtom(new Atom("CA", bs->cacrd()));
		res->addAtom(new Atom("C", bs->ccrd()));
		res->addAtom(new Atom("O", bs->ocrd()));
		resList.push_back(res);
	}

	StructureInfo si(resList);
	si.updateTorsion();
	si.updateSecondaryStructure();
	SasaPSD rsp;
	si.updateSAI(&rsp);
	for(int i=0;i<bsList.size();i++){
		bsList.at(i)->sscode = si.getSS(i);
		bsList.at(i)->data_[0] = si.getPhi(i);
		bsList.at(i)->data_[1] = si.getPsi(i);
		bsList.at(i)->data_[2] = si.getOmg(i);
		bsList.at(i)->data_[3] = si.getSai(i);
	}

	Residue* p2;
	Atom* p3;
	for(int i=0;i<bsList.size();i++){
		p2 = resList.at(i);
		for(int j=0;j<p2->getAtomList()->size();j++)
		{
			p3 = p2->getAtomList()->at(j);
			delete p3;

		}
		delete p2;
	}
}

void proteinchain2BackboneSiteList(ProteinChain* pc, vector<BackBoneSite>& bsList){
	StructureInfo si(pc);
	si.updateTorsion();
	si.updateSecondaryStructure();

	SasaPSD rsp;
	si.updateSAI(&rsp);

//	cout << "pc to bbSitesList" << endl;

	int num = si.getResNum();
	for(int i=0;i<num;i++){
		BackBoneSite bs;
		Residue* res = si.getResidue(i);
	//	cout << "resID: " << res->resID	<< endl;
		bs.pdbid = pc->getPDBID();
		bs.chainid = pc->getChainID();
		bs.resid = si.getResID(i);
		bs.resname = res->getType();
		bs.resseq = i;
		bs.sscode = si.getSS(i);
		bs.data_[0] = si.getPhi(i);
		bs.data_[1] = si.getPsi(i);
		bs.data_[2] = si.getOmg(i);
		bs.data_[3] = si.getSai(i);

	//	cout << "STRUCTUREINFO: PC2BBSITES: " << bs.resid << " " << bs.resname << " " << bs.data_[0] << " " << bs.data_[1] << " " << bs.data_[3] << endl;
		Atom* aN = res->getAtom("N");
		Atom* aCA = res->getAtom("CA");
		Atom* aC = res->getAtom("C");
		Atom* aO = res->getAtom("O");
		XYZ N = aN->getCoord();
		XYZ CA = aCA->getCoord();
		XYZ C = aC->getCoord();
		XYZ O = aO->getCoord();
		bs.data_[4] = N[0];
		bs.data_[5] = N[1];
		bs.data_[6] = N[2];
		bs.data_[7] = CA[0];
		bs.data_[8] = CA[1];
		bs.data_[9] = CA[2];
		bs.data_[10] = C[0];
		bs.data_[11] = C[1];
		bs.data_[12] = C[2];
		bs.data_[13] = O[0];
		bs.data_[14] = O[1];
		bs.data_[15] = O[2];

		bsList.push_back(bs);
	}

}


/*
 * added by xuyang, 2019.3.31
 */

void pdb2BackboneSiteList(PDB* pdb, vector<BackBoneSite*>& bsList, std::string &ssseq) {
    StructureInfo si(pdb);
    si.updateTorsion();
    si.updateSecondaryStructure();

    ssseq.clear();

    SasaPSD rsp;
    si.updateSAI(&rsp);
    int num = si.getResNum();
    for(int i=0;i<num;i++){
        BackBoneSite* bs = new BackBoneSite();
        Residue* res = si.getResidue(i);
        bs->pdbid = pdb->getPDBID();
        bs->chainid = res->chainID;
        bs->resid = si.getResID(i);
        bs->resname = res->getType();
        bs->resseq = i;
        bs->sscode = si.getSS(i);
        bs->data_[0] = si.getPhi(i);
        bs->data_[1] = si.getPsi(i);
        bs->data_[2] = si.getOmg(i);
        bs->data_[3] = si.getSai(i);

        Atom* aN = res->getAtom("N");
        Atom* aCA = res->getAtom("CA");
        Atom* aC = res->getAtom("C");
        Atom* aO = res->getAtom("O");
        XYZ N = aN->getCoord();
        XYZ CA = aCA->getCoord();
        XYZ C = aC->getCoord();
        XYZ O = aO->getCoord();
        bs->data_[4] = N[0];
        bs->data_[5] = N[1];
        bs->data_[6] = N[2];
        bs->data_[7] = CA[0];
        bs->data_[8] = CA[1];
        bs->data_[9] = CA[2];
        bs->data_[10] = C[0];
        bs->data_[11] = C[1];
        bs->data_[12] = C[2];
        bs->data_[13] = O[0];
        bs->data_[14] = O[1];
        bs->data_[15] = O[2];

        bsList.push_back(bs);

        ssseq.push_back(si.getSS(i));
    }
}

} /* namespace NSPdesignseq */
