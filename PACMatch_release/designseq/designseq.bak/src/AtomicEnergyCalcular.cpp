/*
 * AtomicEnergyCalcular.cpp
 *
 *  Created on: 2017Äê12ÔÂ2ÈÕ
 *      Author: notxp
 */

#include "designseq/AtomicEnergyCalcular.h"

namespace NSPdesignseq {

AtomicEnergyCalcular::AtomicEnergyCalcular() {
	// TODO Auto-generated constructor stub
	this->lamdaRep = 0.1;
	this->lamdaAtr = 0.15;
	this->wdNP_NP = -1.3;
	this->wdP_NP = -0.7;
	this->wdP_P = 0.0;
	this->wdAro = -1.5;

	this->wdHB = -0.2;
	this->lamdaHBRep = 0.1;
	this->lamdaHBAtr = 0.15;

	initVdwEnergyTable();
}

AtomicEnergyCalcular::AtomicEnergyCalcular(DesignParameters* para){
	this->lamdaRep = para->vdwLamda1;
	this->lamdaAtr = para->vdwLamda2;
	this->wdNP_NP = para->wdNP_NP;
	this->wdP_NP = para->wdP_NP;
	this->wdP_P = para->wdP_P;
	this->wdHB = para->wdHB;
	this->wdAro = para->aroWD;

	this->lamdaHBRep = para->vdwLamda1;
	this->lamdaHBAtr = para->vdwLamda2;
	this->hbAngleFactor = para->hbAngleFactor;

	this->para = para;

	if(abs(lamdaRep) < 0.001)
		lamdaRep = 0.001;
	if(abs(lamdaAtr) < 0.001)
		lamdaAtr = 0.001;
	if(abs(lamdaHBRep) < 0.001)
		lamdaHBRep = 0.001;
	if(abs(lamdaHBAtr) < 0.001)
		lamdaHBAtr = 0.001;

	initVdwEnergyTable();
}

void AtomicEnergyCalcular::initVdwEnergyTable()
{
	int i,j;
	float dd,d,d0,u;
	for(i=0;i<100;i++){
		d0 = 3.0+i*0.01 + 0.005;

		for(j=0;j<3200;j++){
			dd = j*0.005 + 0.0025;
			d = sqrt(dd);
			u = (d-d0)/d0/lamdaRep;
			vdwEnergyTable1[i][j] = u*u;
		}
	}

	for(i=0;i<100;i++){
		d0 = 3.0+i*0.01 + 0.005;
		for(j=0;j<10000;j++){
			dd = j*0.005 + 0.0025;
			d = sqrt(dd);
			u = (d-d0)/d0/lamdaAtr;
			vdwEnergyTable2[i][j] = exp(-u*u);
		}
	}

}


/*
float AtomicEnergyCalcular::getEnergy(Tern* coordA, Tern* dirA, AtomProperty* apA,  Tern* coordB, Tern* dirB, AtomProperty* apB)
{
	if(!coordA->neighborTo(*coordB,6.0))
		return 0.0;
	float wd = 0, rescale;
	bool hbonded = false;
	int neighbor = apA->connectNum + apB->connectNum -1;
	if(neighbor  == 1)
		rescale = 1.5;
	else
		rescale = neighbor;

	if(!apA->isPolar && !apB->isPolar)
		wd = wdNP_NP;
	else if((apA->isHDonor && apB->isHAcceptor) || (apA->isHAcceptor && apB->isHDonor))
	{
		wd = wdHB;
		hbonded = true;
		return 0.0;
	}
	else if(apA->isPolar && apB->isPolar)
		wd = wdP_P;
	else if(apA->isPolar || apB->isPolar)
		wd = wdP_NP;

	Tern ab = *coordB - *coordA;
	float d = ab.length();

	float r1, r2, ang1, ang2, cosAngle;
	r1 = apA->vdwRadius;
	r2 = apB->vdwRadius;

	if(apA->isAromatic)
	{

		if(abs(dirA->length() - 1) > 0.0001)
			cerr << "invalid dir: " << dirA->t[0] << " " << dirA->t[1] << " " << dirA->t[2] << endl;
		cosAngle = ab*(*dirA)/d;
		if(cosAngle > 1)
			cosAngle = 1;
		if(cosAngle < -1)
			cosAngle = -1;
		if(cosAngle > 1 || cosAngle < -1)
		{
			cerr << "invalid cosAngle: " << cosAngle << endl;
			cerr << "dir: " << dirA->x() << " " << dirA->y() << " " << dirA->z() << endl;
		}
		ang1 = acos(cosAngle)*180/PI;


		r1 = apA->getSp2Radius(ang1);
	//	cout << "Ar1: " << r1 << " Ang1: " << ang1 << endl;
	}

	if(apB->isAromatic)
	{
		if(abs(dirB->length() - 1) > 0.01)
			cerr << "invalid dir: " << dirB->t[0] << " " << dirB->t[1] << " " << dirB->t[2] << endl;
		cosAngle = ab*(*dirB)/d;
		if(cosAngle > 1)
			cosAngle = 1;
		if(cosAngle < -1)
			cosAngle = -1;
		if(cosAngle > 1 || cosAngle < -1)
		{
			cerr << "invalid cosAngle: " << cosAngle << endl;
			cerr << "ab: " << ab.x() << " " << ab.y() << " " << ab.z() << endl;
			cerr << "dir: " << dirB->x() << " " << dirB->y() << " " << dirB->z() << endl;
		}
		ang2 = acos(cosAngle)*180/PI;
		r2 = apB->getSp2Radius(ang2);

	}

	if(hbonded)
	{
		r1 = apA->polarRadius;
		r2 = apB->polarRadius;
		return getEnergy(d,r1+r2,wd,1.0);
	}
	else
		return getEnergy(d,r1+r2,wd,rescale);
}
*/

float AtomicEnergyCalcular::getAtomEnergy(XYZ* coordA, XYZ* dirA, AtomProperty* apA,  XYZ* coordB, XYZ* dirB, AtomProperty* apB)
{
	float dd = coordA->squaredDistance(*coordB);
	float d = sqrt(dd);
	if((apA->isHDonor && apB->isHAcceptor) || (apA->isHAcceptor && apB->isHDonor))
		return getHBEnergy(d, 2.8, this->wdHB);

	if(dd > 49)
		return 0.0;

	float wd = 0, rescale;
	int neighbor = apA->connectNum + apB->connectNum -1;
	if(neighbor == 1)
		rescale = 1.0;
	else
		rescale = neighbor/1.5;


	if(!apA->isPolar && !apB->isPolar)
		wd = wdNP_NP;
	else if(apA->isPolar && apB->isPolar)
		wd = wdP_P;
	else if(apA->isAromatic || apB->isAromatic)
		wd = wdAro;
	else
		wd = wdP_NP;

	XYZ ab = *coordB - *coordA;

	float r1, r2, cosAngleSquare, tmp;
	r1 = apA->vdwRadius;
	r2 = apB->vdwRadius;


	if(apA->isAromatic)
	{
		tmp = ab * (*dirA);
		cosAngleSquare = tmp*tmp/dd;
		r1 = apA->getSp2Radius2(cosAngleSquare);
	}

	if(apB->isAromatic)
	{
		tmp = ab*(*dirB);
		cosAngleSquare = tmp*tmp/dd;
		r2 = apB->getSp2Radius2(cosAngleSquare);
	}


	return getEnergy2(dd,r1+r2,wd,rescale);
}


AtomicEnergyCalcular::~AtomicEnergyCalcular() {
	// TODO Auto-generated destructor stub
}

} /* namespace xenergy */

