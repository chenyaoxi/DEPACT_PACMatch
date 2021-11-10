/*
 * testrotation.cpp
 *
 *  Created on: 2016年9月29日
 *      Author: hyliu
 */


#include "geometry/rigidtransform.h"
using namespace NSPgeometry;
int main(int argc,char **argv){
	Rotation r1(QuaternionCrd(XYZ(1,2,3),30.0),XYZ());
	RigidTransform rt(QuaternionCrd(XYZ(1,2,3),30.0),XYZ(),XYZ(3,-1,4));
	XYZ pa(3.2,4.3,1.5);
	XYZ pbb=r1.applytoCopy(pa);
	XYZ pb=rt.applytoCopy(pa);
	std::cout <<pb.toString() <<std::endl;
	std::cout <<pbb.toString() <<std::endl;
	Rotation r2(QuaternionCrd(XYZ(1,0,-3),60.0),XYZ(1,3,5));
	RigidTransform rtnew=applyRotation(r2,rt);
	XYZ pc=r2.applytoCopy(pb);
	XYZ pcc=rtnew.applytoCopy(pa);
	std::cout <<pc.toString() <<std::endl;
	std::cout <<pcc.toString() <<std::endl;


/*	XYZ p1(3.0,-5.4,8.3);
	XYZ p2(4.0,-7.4,9.3);
	XYZ f1(-2,3,0.5);
	XYZ f2(0,4.1,-0.5);
	Rotation R=rotationaligntwovectors(p1,p2-p1,f2-f1);
	XYZ newf2=p1+f2-f1;
	R.apply(&newf2);
	std::cout <<"p2: " <<(p2-p1).toString() <<std::endl;
	std::cout <<"newf2: " <<(newf2-p1).toString() <<std::endl;

	XYZ r1(0.2,0.3,0.5);
	double theta1=sqrt(r1.squarednorm())/3.14159265*180.0;
	QuaternionCrd q1(r1,theta1);
	XYZ r2(0.1,-0.2,0.4);
	double theta2=sqrt(r2.squarednorm())/3.14159265*180.0;
	QuaternionCrd q2(r2,theta2);
	QuaternionCrd q3=q2*q1;
	std::cout <<q3.diff(q2) <<"   " <<theta1<<std::endl;
	Rotation m;
	m.init(q2,XYZ(0,0,0));
	std::cout <<r1.toString() <<std::endl;
	m.apply(&r1);
	std::cout <<r1.toString() <<std::endl;
	XYZ rt=r1+r2;
	double theta4=sqrt(rt.squarednorm())/3.14159265*180.0;
	QuaternionCrd q4(rt,theta4);
	std::cout <<q3.diff(q4)<<std::endl;
	*/
}

