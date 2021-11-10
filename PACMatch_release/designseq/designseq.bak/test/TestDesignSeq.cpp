/*
 * TestDesignSeq.cpp
 *
 *  Created on: 2017��12��27��
 *      Author: notxp
 */


#include <iostream>
#include <vector>
#include <string>
#include <time.h>
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
#include "designseq/SeqMinimize.h"
#include "designseq/CmdArgs.h"

using namespace std;
using namespace NSPdesignseq;
using namespace NSPgeometry;

int main(int argc, char** args){
	/*
	 * input is PDB structure
	 */

	CmdArgs cmd(argc, args);

	if(!cmd.specifiedOption("-in") || !cmd.specifiedOption("-out")){
		cout << "Usage:\nDesignSeq -in $INPUTPDB -out $OUTPUTPDB" << endl;
		cout << "DesignSeq -in $INPUTPDB -out $OUTPUTPDB -n $SEQNUM" << endl;
		cout << "DesignSeq -in $INPUTPDB -out $OUTPUTPDB -resFile $RESFILE" << endl;
	}

	int n = 1;
	if(cmd.specifiedOption("-n")){
		n = atoi(cmd.getValue("-n").c_str());
	}



	clock_t start = clock();


	string s = cmd.getValue("-in");
	string pdbID = "unk";
	string output = cmd.getValue("-out");

	if(output.length() < 4 || output.substr(output.length()-4, 4) != ".pdb")
		output = output + ".pdb";


	PDB pdb(s, pdbID);
	ProteinChain* pc = pdb.getFirstChain();
	string natSeq = pc->getSequence();

	vector<BackBoneSite> bsList0;
	proteinchain2BackboneSiteList(pc, bsList0);


	vector<BackBoneSite*> bsList;
	for(int i=0;i<bsList0.size();i++){
		bsList.push_back(&bsList0.at(i));
	}


	DesignParameters dp;
	if(cmd.specifiedOption("-paraFile")){
		std::string parafile=cmd.getValue("-paraFile");
		dp=DesignParameters(parafile);
	}
	S1EnergyTable s1Etable;
	S2MatrixFinder s2Etable;


	DesignTemplate* dt;

	if(cmd.specifiedOption("-resFile")){
		string resFile = cmd.getValue("-resFile");
		dt = new DesignTemplate(bsList, resFile, &dp, s1Etable, s2Etable);
	}
	else{
		dt = new DesignTemplate(bsList, &dp, s1Etable, s2Etable);
	}


//	dt.checkResInvolvedMap();
	dt->loadS1S2(s1Etable, s2Etable);
//	dt.printRotChoice();
	dt->loadSingleResidueEnergy();
	dt->loadPairwiseEnergy();

	DesignMC mc(dt);
	RotSequence unit(dt->resNum);


	for(int i=1;i<=n;i++)
	{
		char x[10];
		sprintf(x,"000%d",i);
		string xx(x);
		xx = xx.substr(xx.length()-3,3);
		string fileName = output;
		if(n > 1)
			fileName = fileName.substr(0, output.length()-4) + "-" + xx + ".pdb";
		mc.mcRun(&unit);
		mc.printPDB(&unit, fileName);
	}


	delete dt;

	clock_t end = clock();
	cout << "time: " << (float)(end-start)/CLOCKS_PER_SEC << endl;

}

