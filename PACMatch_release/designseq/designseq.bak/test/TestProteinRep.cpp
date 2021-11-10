/*
 * TextProteinRep.cpp
 *
 *  Created on: 2017Äê10ÔÂ31ÈÕ
 *      Author: notxp
 */

#include <iostream>
#include <vector>
#include <string>
#include "designseq/ProteinRep.h"
#include "designseq/StructureInfo.h"
#include "designseq/S1EnergyTable.h"
#include "geometry/xyz.h"
#include "designseq/AtomLib.h"
#include "designseq/StructureInfo.h"
#include "designseq/S2MatrixFinder.h"
#include "designseq/AtomicEnergyCalcular.h"


using namespace std;
using namespace NSPdesignseq;
using namespace NSPgeometry;

int main(int argc, char** args){

	string s(args[1]);
	string pdbID = "unk";
	PDB pdb(s, pdbID);
	ProteinChain* pc = pdb.getFirstChain();
	Residue* res = pc->getResidue("10");
	vector<Residue*> resList = pc->getResList();
	vector<XYZ> caList1;
	for(Residue* resx : resList){
		caList1.push_back(resx->getAtom("CA")->getCoord());
	}

	string s2(args[2]);
	PDB pdb2(s2,pdbID);
	ProteinChain* pc2 = pdb2.getFirstChain();
	vector<Residue*> res2List = pc2->getResList();
	vector<XYZ> caList2;
	for(Residue* resx : res2List){
		caList2.push_back(resx->getAtom("CA")->getCoord());
	}

	vector<double> weights;
	for(int i=0;i<res2List.size();i++){
		weights.push_back(1.0);
	}
	QuatFit qf;
	double rmsd1ubq = qf.setup(caList1, caList2, weights);
    cout << "1ubq RMSD : " << rmsd1ubq << endl;


	LocalFrame cs = res->getCoordSystem();
	XYZ x = cs.axis_[0];
	XYZ y = cs.axis_[1];
	XYZ z = cs.axis_[2];


	printf("%8.3f %8.3f %8.3f\n",x[0],y[0],z[0]);
	printf("%8.3f %8.3f %8.3f\n",x[1],y[1],z[1]);
	printf("%8.3f %8.3f %8.3f\n",x[2],y[2],z[2]);

	printf("%8.3f %8.3f %8.3f\n",cs.origin_[0], cs.origin_[1], cs.origin_[2]);

	AtomLib atLib;
	int id = atLib.uniqueNameToID("ALA-CB");
	printf("ALA-CB %3d\n",id);

	cout << "test structure info:" << endl;

	StructureInfo si(&pdb);
	SasaPSD sp;


	si.updateTorsion();
	si.updateSecondaryStructure();
	si.updateSAI(&sp);


	int num = si.getResNum();
	for(int i=0;i<num;i++){
		float sai = si.getSai(i);
		cout << i << " " << sai << endl;
	}

	cout << "start " << endl;

	S1EnergyTable s1ET;
	for(float sai=0;sai < 1; sai += 0.02){
		int id = s1ET.saiToInt(sai);
		cout << sai << " " << id << endl;
	}

	vector<BackBoneSite> bbSites;
	cout << bbSites.size() << endl;

	proteinchain2BackboneSiteList(pc, bbSites);

/*

	for(int i=0;i<bbSites.size();i++)
	{

		BackBoneSite& bs = bbSites.at(i);

		printf("%-2d %s ", bs.resid, bs.resname.c_str());
		//cout << "start search s1:" << endl;
		float* s1 = s1ET.getS1(bs);
		//cout << "search s1 finished" << endl;
		for(int k=0;k<20;k++){

			printf("%6.4f ", s1[k]);
		}
		cout << endl;
	}


*/

	Residue* res1 = pc->getResidue("5");
	cout << "res1: " << endl;
	Residue* res2 = pc->getResidue("30");
	cout << "res2: " << endl;

	LocalFrame cs1 = res1->getCoordSystem();
	cout << "cs1" << endl;
	printf("ori: %7.3f %7.3f %7.3f\n", cs1.origin_[0], cs1.origin_[1], cs1.origin_[2]);
	printf("  X: %7.3f %7.3f %7.3f\n", cs1.axis_[0][0], cs1.axis_[1][0], cs1.axis_[2][0]);
	printf("  Y: %7.3f %7.3f %7.3f\n", cs1.axis_[0][1], cs1.axis_[1][1], cs1.axis_[2][1]);
	printf("  Z: %7.3f %7.3f %7.3f\n", cs1.axis_[0][2], cs1.axis_[1][2], cs1.axis_[2][2]);
	cout << endl;


	LocalFrame cs2 = res2->getCoordSystem();
	cout << "cs2" << endl;

	Residue* res3 = pc->getResidue("27");

	cout << "res3" << endl;
	Residue* res4 = pc->getResidue("43");
	cout << "res4" << endl;
	LocalFrame cs3 = res3->getCoordSystem();
	cout << "cs3" << endl;
	LocalFrame cs4 = res4->getCoordSystem();
	cout << "cs4" << endl;

	ResPairOrientation rp1(cs1, cs2);
	cout << "rp1" << endl;
	ResPairOrientation rp2(cs3, cs4);
	cout << "rp2" << endl;

	cout << rp1.toString() << endl;
	cout << rp2.toString() << endl;

	double rms = rp1.rmsd(rp2);




	printf("rmsd2: %6.4f\n", rms);

	cout << "test S2 matrix finder" << endl;

	S2MatrixFinder s2et;
	cout << "s2 etable inited " << endl;


	BackBoneSite bsA = bbSites.at(4);
	BackBoneSite bsB = bbSites.at(24);

	BackboneSitesPair rp(&bsA, &bsB);
	AAScoreMatrix sm;

	cout << "pair info" << endl;
	cout << rp.siteA->sscode << endl;
	cout << rp.siteA->data_[3] << endl;

	cout << rp.siteB->sscode << endl;
	cout << rp.siteB->data_[3] << endl;

	cout << rp.ori.toString() << endl;

	s2et.getSM(&rp, &sm);

	for(int i=0;i<20;i++){
		for(int j=0;j<20;j++){
			double s = sm.getValue(i,j);
			printf("%6.3f ",s);
		}
		cout << endl;
	}



}


