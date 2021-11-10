/*
 * AtomicEnergyCalcular.h
 *
 *  Created on: 2017Äê12ÔÂ2ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_ATOMICENERGYCALCULAR_H_
#define DESIGNSEQ_ATOMICENERGYCALCULAR_H_


#include <iostream>
#include <map>
#include "designseq/DesignParameters.h"
#include "designseq/ProteinRep.h"
#include "geometry/xyz.h"

namespace NSPdesignseq {

using namespace std;



class AtomicEnergyCalcular {
private:

	float lamdaRep;
	float lamdaAtr;
	float wdNP_NP;
	float wdP_NP;
	float wdP_P;
	float wdAro;


	float wdHB;

	float lamdaHBRep;
	float lamdaHBAtr;
	float hbAngleFactor;

	/*
	 * non-hbond
	 * first dimension: d0, from 3.0 to 4.0
	 * second dimension: d^2, from 0 to 50
	 */
	float vdwEnergyTable1[100][3200];
	float vdwEnergyTable2[100][10000];

	/*
	 * hbond
	 * first dimension: d0, from 2.5 to 3.5
	 * second dimension: d^2, from 0 to 30
	 */

	//float hbEnergyTable[100][3000];


public:

	DesignParameters* para;
	AtomicEnergyCalcular();
	AtomicEnergyCalcular(DesignParameters* para);
	void initVdwEnergyTable();
	float getEnergy(float d, float d0, float wd, float rescale){



		float e;
		if(d < d0)
		{
			float u = (d-d0)/d0/lamdaRep;
			return u*u+wd/rescale;
		}
		else
		{
			float u = (d-d0)/d0/lamdaAtr;
			if(u>4)
				e = 0.0;
			else
				e = wd*exp(-u*u);
			return e/rescale;
		}
	}

	float getEnergy2(float dd, float d0, float wd, float rescale){


		float d0Square = d0*d0;

		int ddIndex = (int)(dd*200);
		if(ddIndex >= 10000)
		{
			cerr << "invalid dd: " << dd << endl;
			exit(1);
		}
		int d0Index = (int)(d0*100-300);
		if(d0Index >= 100)
		{
			cerr << "invalid d0: " << d0 << endl;
			exit(1);
		}
		if(dd < d0Square)
		{
			return vdwEnergyTable1[d0Index][ddIndex] + wd/rescale;
		}
		else
		{
			return wd/rescale * vdwEnergyTable2[d0Index][ddIndex];
		}
	}

	float getHBEnergy(float d, float d0, double wd){
			float e;
			wd = wd*1.5;
			if(d < d0)
			{
				float u = (d-d0)/d0/lamdaHBRep;
				return u*u+wd;
			}
			else
			{
				float u = (d-d0)/d0/lamdaHBAtr;
				if(u>4)
					e = 0.0;
				else
					e = wd*expf(-u*u);
				return e;
			}
	}


	float getAtomEnergy(XYZ* coordA, XYZ* dirA, AtomProperty* apA,  XYZ* coordB, XYZ* dirB, AtomProperty* apB);



	virtual ~AtomicEnergyCalcular();
};



} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_ATOMICENERGYCALCULAR_H_ */
