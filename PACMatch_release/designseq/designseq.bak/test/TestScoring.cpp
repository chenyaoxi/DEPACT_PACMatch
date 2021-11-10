/*
 * TestScoring.cpp
 *
 *  Created on: 2017��12��25��
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
#include "designseq/DesignTemplate.h"
#include "designseq/AAProbabilityArray.h"

using namespace std;
using namespace NSPdesignseq;
using namespace NSPgeometry;

int main(int argc, char** args){

	/*
	 * input is PDB structure
	 */
	if(argc != 2)
	{
		cout << "testScoring $PDBFile" << endl;
	}
	string s(args[1]);
	string pdbID = "unk";
	PDB pdb(s, pdbID);
	ProteinChain* pc = pdb.getFirstChain();

	/*
	 * calculate backbone torsion angle, secondary structure, and SASA
	 */
	StructureInfo strInfo(pc);

	SasaPSD sasaPSD;
	strInfo.updateTorsion();
	strInfo.updateSecondaryStructure();
	strInfo.updateSAI(&sasaPSD);



	/*
	 * print local structure information
	 */
	cout << "local structure information:" << endl;
	int resNum = pc->getChainLength();
	for(int i=0;i<resNum;i++){
		Residue* res = pc->getResList().at(i);
		string resID = res->resID;
		char ss = strInfo.getSS(i);
		float phi = strInfo.getPhi(i);
		float psi = strInfo.getPsi(i);
		float omg = strInfo.getOmg(i);
		float sai = strInfo.getSai(i);
		string triName = res->triName;
		printf("%3s %s %c %7.2f %7.2f %7.2f %5.2f\n", resID.c_str(), triName.c_str(), ss, phi, psi, omg, sai);
	}


	/*
	 * convert ProteinChain to BackboneSite list
	 * in this step, the torsion angles, secondary structure, and SASA will be calculated
	 */

	vector<BackBoneSite> bsList;
	proteinchain2BackboneSiteList(pc, bsList);


	/*
	 * calculate S1 score
	 */

	S1EnergyTable s1Etable;
	ResName rn;
	cout << "S1 score: " << endl;

	for(int i=0;i<resNum;i++){
		Residue* res = pc->getResList().at(i);
		string resID = res->resID;
		string triName = res->triName;

		AAProbabilityArray pa;
		s1Etable.getS1(bsList.at(i), &pa);
		int intName = rn.triToInt(triName);
		float s1 = pa.getScore(intName);

		printf("%3s %s %6.3f\n", resID.c_str(), triName.c_str(), s1);
	}

	cout << "S2 score: " << endl;


	/*
	 * update residue pairs
	 */
	S2MatrixFinder s2Etable;
	for(int i=0;i<resNum;i++){
		BackBoneSite* bsA = &(bsList.at(i));
		string triNameA = pc->getResList().at(i)->triName;
		int intNameA = rn.triToInt(triNameA);
		string resIDA = pc->getResList().at(i)->resID;

		for(int j=i+1;j<resNum;j++){
			BackBoneSite* bsB = &(bsList.at(j));
			string triNameB = pc->getResList().at(j)->triName;
			int intNameB = rn.triToInt(triNameB);
			string resIDB = pc->getResList().at(j)->resID;
			float cbDistance = bsA->cbcrd().distance(bsB->cbcrd());
			/*
			 * residue pair is defined by CB atom distance smaller than 8.0
			 */
			if(cbDistance > 8.0) continue;
			BackboneSitesPair pair(&(bsList.at(i)), &(bsList.at(j)));
			AAScoreMatrix sm;
			s2Etable.getSM(&pair, &sm);

			float s2 = sm.getValue(intNameA, intNameB);
			printf("%3s %s %3s %s %7.3f\n",resIDA.c_str(), triNameA.c_str(), resIDB.c_str(), triNameB.c_str(), s2);
		}
	}



	/*
	 * calculate atomic score
	 */

	cout << "atomic energy:" << endl;

	AtomicEnergyCalcular ec;

	AtomLib atLib;
	vector<Rotamer> rotList;
	for(int i=0;i<resNum;i++){
		Residue* res = pc->getResList().at(i);
		Rotamer rot = res->natRotamer(&atLib);
		rotList.push_back(res->natRotamer(&atLib));
	}

	for(int i=0;i<resNum;i++){
		string resIDA = pc->getResList().at(i)->resID;
		BackBoneSite* bsA = &(bsList.at(i));
		Rotamer* rotA = &(rotList.at(i));
		string triNameA = rotA->triName;

		Conformer* confA = new Conformer(rotA, bsA, &atLib);


		for(int j=i+2;j<resNum;j++){
			string resIDB = pc->getResList().at(j)->resID;


			BackBoneSite* bsB = &(bsList.at(j));
			Rotamer* rotB = &(rotList.at(j));
			string triNameB = rotB->triName;
			Conformer* confB = new Conformer(rotB, bsB, &atLib);

			int seqSep = j-i;
			if(seqSep > 5) seqSep = 5;
		//	cout << "calculate pair energy: " << endl;
			float e = pairEnergy(confA, confB, seqSep, 1.0, &ec);

			if(e != 0){
				printf("%3s %3s %3s %3s %7.3f\n", resIDA.c_str(), triNameA.c_str(), resIDB.c_str(), triNameB.c_str(), e);
			}
		}
	}

}


