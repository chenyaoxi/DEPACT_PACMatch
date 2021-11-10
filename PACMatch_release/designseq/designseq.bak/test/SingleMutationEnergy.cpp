/*
 * SingleMutationEnergy.cpp
 *
 *  Created on: 2018Äê2ÔÂ26ÈÕ
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
	if(argc != 5)
	{
		cout << "singleMut $PDBFile $ChainID $ResID $MutType" << endl;
	}
	string s(args[1]);
	string pdbID = "unk";
	PDB pdb(s, pdbID);

	char chainID = args[2][0];
	string resid(args[3]);
	string mutType(args[4]);

	ProteinChain* pc = pdb.getChain(chainID);
	if(pc == NULL){
		cerr << "invalid chainID: " << chainID << endl;
		exit(0);
	}

	vector<BackBoneSite> bsList0;
	proteinchain2BackboneSiteList(pc, bsList0);


	vector<BackBoneSite*> bsList;
	for(int i=0;i<bsList0.size();i++){
		bsList.push_back(&bsList0.at(i));
	}



	Residue* res = pc->getResidue(resid);
	if(res == NULL){
		cerr << "invalid resid: " << resid << endl;
		exit(0);
	}

	int resSeqID = res->resSeqID;






}

