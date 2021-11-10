/*
 * SeqMinimize.cpp
 *
 *  Created on: 2017Äê12ÔÂ26ÈÕ
 *      Author: notxp
 */

#include "designseq/SeqMinimize.h"

namespace NSPdesignseq {

	RotSequence::RotSequence(int len){
		this->seqLength = len;
		this->seqChoices = new int[len];
		for(int i=0;i<len;i++)
			this->seqChoices[i] = 0;
	}

	RotSequence::RotSequence(vector<RotamerGroup*>& groups){
		int len = groups.size();
		this->seqLength = len;
		this->seqChoices = new int[len];

		for(int i=0;i<len;i++){
			this->rotGroups.push_back(groups.at(i));
		}
		setRandomChoice();
	}

	void RotSequence::copyValue(RotSequence* from){
		if(this->seqLength != from->seqLength){
			cerr << "rotamer sequence length not equal" << endl;
			exit(1);
		}

		this->rotGroups.clear();
		for(int i=0;i<seqLength;i++){
			this->seqChoices[i] = from->seqChoices[i];
			this->rotGroups.push_back(from->rotGroups.at(i));
		}
	}

	void RotSequence::setRandomChoice(){
	//	cout << "set Random choice:" << endl;
		for(int i=0;i<this->seqLength;i++){
			int choiceNum = this->rotGroups.at(i)->rotNum;
	//		cout << "rot number: " << choiceNum << endl;
			int randNum = rand()%choiceNum;
			this->seqChoices[i] = randNum;
		}
	}

	void RotSequence::applyMutation(int pos, int choice){
		this->seqChoices[pos] = choice;
	}

	string RotSequence::toAASequence(){
		ResName rn;
		char seq[this->seqLength+1];
		for(int i=0;i<seqLength;i++){
			Rotamer* rot = this->rotGroups[i]->rotList.at(this->seqChoices[i]);
			seq[i] = rn.triToSin(rot->triName);
		}
		seq[this->seqLength] = '\0';
		string s = seq;
		return s;
	}

	double RotSequence::mutEnergy(int pos, int mutChoice, EnergyCalculatorTemplate* ec){
		int oldChoice = this->seqChoices[pos];
		return resInvolvedEnergy(pos, mutChoice, ec) - resInvolvedEnergy(pos, oldChoice, ec);
	}

	double RotSequence::resInvolvedEnergy(int pos, int choice, EnergyCalculatorTemplate* ec){
		double e1 = ec->eaList.at(pos)->getEnergy(choice);
		double e2 = 0;
		int pairNum = ec->involvedPairs.at(pos)->size();
		for(int i=0;i<pairNum;i++){
			int pairID = ec->involvedPairs.at(pos)->at(i);
			BackboneSitesPair* pair = ec->bsPairs.at(pairID);
			EnergyMatrix* em = ec->emList.at(pairID);
			int posA = pair->siteA->resseq;
			int posB = pair->siteB->resseq;
			if(pos == posA){
				e2 += em->getEnergy(choice,this->seqChoices[posB]);
			}
			else if(pos == posB){
				e2 += em->getEnergy(this->seqChoices[posA], choice);
			}
			else{
				cerr << "posA: " << posA << " posB: " << posB << " mutPos: " << pos << " resInvolved Energy error" << endl;
				exit(1);
			}
		}
		return e1+e2;
	}

	void RotSequence::acceptMut(int pos, int mutChoice){
		this->seqChoices[pos] = mutChoice;
	}

	double RotSequence::totEnergy(EnergyCalculatorTemplate* ec){
		double e = 0;
		for(int i=0;i<ec->eaList.size();i++){
			//cout << "cal totEnergy: " << i << endl;
			e += ec->eaList.at(i)->getEnergy(this->seqChoices[i]);
		}
		for(int i=0;i<ec->emList.size();i++){
			BackboneSitesPair* pair = ec->bsPairs.at(i);
			EnergyMatrix* em = ec->emList.at(i);
			int posA = pair->siteA->resseq;
			int posB = pair->siteB->resseq;
			e += em->getEnergy(this->seqChoices[posA], this->seqChoices[posB]);
		}
		return e;
	}

	void RotSequence::printDetailEnergy(EnergyCalculatorTemplate* ec){
		double e = 0;
		double e1 = 0;
		double e2 = 0;
		for(int i=0;i<ec->eaList.size();i++){
			Rotamer* rot = getRotamer(i);
			e = ec->eaList.at(i)->getEnergy(getChoice(i));
			e1 += e;
			printf("Res: %3d Rot: %7s Energy: %10.3f\n", i+1, rot->rotName.c_str(), e);
		}

		for(int i=0;i<ec->emList.size();i++){
			BackboneSitesPair* pair = ec->bsPairs.at(i);
			EnergyMatrix* em = ec->emList.at(i);
			int posA = pair->siteA->resseq;
			int posB = pair->siteB->resseq;
			Rotamer* rotA = getRotamer(posA);
			Rotamer* rotB = getRotamer(posB);
			e = em->getEnergy(this->seqChoices[posA], this->seqChoices[posB]);
			e2+= e;
			printf("Pair: %3d %3d Rot: %-7s %-7s Energy: %10.3f\n", posA+1, posB+1, rotA->rotName.c_str(), rotB->rotName.c_str(), e );
		}

		printf("TotalEnergy: E1: %7.3f E2: %7.3f\n", e1, e2);
	}

	RotSequence::~RotSequence(){
		delete[] this->seqChoices;
	}

    bool DesignMC::accept(double mutEnergy, double T){
    	if(T == 9){
    		if(mutEnergy > 0)
    			return false;
    		else
    			return true;
    	}
    	else if(mutEnergy > 0){
    		float pAc = exp(-1.0*mutEnergy/T);
    		float r = 1.0*rand()/RAND_MAX;
    		return r < pAc;
    	}
    	else
    		return true;
    }

    void DesignMC::mcRun(RotSequence* result){
         RotSequence unit(dt->rotGroups);
         result->copyValue(&unit);


     //    cout << "start MC: " << endl;
        // unit.checkChoice();

     //    unit.printDetailEnergy(dt);

         double energy = unit.totEnergy(dt);
    //     cout << "totEnergy: " << energy << endl;
         int len = dt->resNum;

         int randPos, randChoice;
         double mutEnergy;

         srand((unsigned)time(NULL));
         double minEnergy = energy;
         int count = 0;
         for(double T = T0;T>T1;T=T*annealFactor){
        	 bool noChange = true;
        	 int acNum = 0;
        	 for(int i=0;i<step;i++){
        		 count ++;
        	//	 if(count %1000 == 0){
        	//		 printf("%-10d %10.4f\n",count,energy);
        	//	 }
        		 randPos = rand()%len;
        		 randChoice = rand()%(unit.choiceNum(randPos));
        		// cout << "pos: " << randPos << " randChoice: " << randChoice << " choiceNum: " << dt->eaList.at(randPos)->getChoiceNum() << endl;
        		 mutEnergy = unit.mutEnergy(randPos,randChoice, dt);
        		 if(mutEnergy < 0)
        			 noChange = false;
        		 if(accept(mutEnergy, T)){
        			 acNum++;
        			 unit.acceptMut(randPos, randChoice);
        			 energy += mutEnergy;
        			 if(energy < minEnergy){
        				 minEnergy = energy;
        				 result->copyValue(&unit);
        			 }
        		 }
        		// cout << "step = " << i << endl;
        	 }
        	// printf("%7.4f %10.4f %10.4f %4d\n",T,energy,minEnergy,acNum);
        	 if(noChange) break;
         }
    }

    void DesignMC::printPDB(RotSequence* unit, string outputFile){
    	ProteinChain pc;
    	char s[20];
    	for(int i=0;i<dt->bsList.size();i++){
    		BackBoneSite* bs = dt->bsList.at(i);
    		Rotamer* rot = unit->getRotamer(i);
    		sprintf(s,"%d",bs->resid);
    		Residue* res = new Residue(string(s), bs->chainid, rot->triName);
    		res->addAtom(new Atom("N", bs->ncrd()));
    		res->addAtom(new Atom("CA", bs->cacrd()));
    		res->addAtom(new Atom("C", bs->ccrd()));
    		res->addAtom(new Atom("O", bs->ocrd()));
    		res->buildRotamer(rot);
    		pc.addResidue(res);
    	}
    	ofstream output(outputFile, ios::out);
    	if(!output.is_open())
    	{
    		cout << "fail to open file " << outputFile << endl;
    		exit(1);
    	}
    	pc.printPDBFormat(output, 1);

    	Residue* p1;
    	Atom* p2;
    	for(int i=0;i<pc.getChainLength();i++){
    		p1 = pc.getResList().at(i);
    		for(int j=0;j<p1->getAtomList()->size();j++){
    			p2 = p1->getAtomList()->at(j);
    			delete p2;
    		}
    		delete p1;
    	}
    	output.close();
    }

    DesignMC::~DesignMC(){

    }



} /* namespace NSPdesignseq */
