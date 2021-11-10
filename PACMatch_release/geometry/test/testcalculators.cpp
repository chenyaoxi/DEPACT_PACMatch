/*
 * testcalculators.cpp
 *
 *  Created on: 2017年10月29日
 *      Author: hyliu
 */

#include "geometry/calculators.h"
#include <iostream>
using namespace NSPgeometry;

int main(){
	std::vector<XYZ> atoms{{1.2,0.8,0.8},{1.2,0.8,0.4},{1.0,0.3,0.1},{0.0,0.0,0.0}};
	std::vector<XYZ> derivb,deriva,derivt;
	std::cout << "bond 1-2: "<< distance(atoms[0],atoms[1],&derivb)<<std::endl;
	std::vector<XYZ> dr={{0.0001,0.0,0.0},{0.0,0.0001,0.0},{0.0,0.0,0.0001}};
	for(int i=0;i<2;++i){
		for(int m=0; m<3;++m) {
			atoms[i] =atoms[i]+dr[m];
			double dp=distance(atoms[0],atoms[1]);
			atoms[i] = atoms[i]-2.0*dr[m];
			double dm=distance(atoms[0],atoms[1]);
			std::cout <<i<<"\t"<<m<<"\t"<<derivb[i][m]<<"\t"<<(dp-dm)/0.0002<<std::endl;
			atoms[i] =atoms[i]+dr[m];
		}
	}
	std::cout << "ang 1-2-3: "<< cos_angle(atoms[0],atoms[1],atoms[2],&deriva)<<std::endl;
	for(int i=0;i<3;++i){
			for(int m=0; m<3;++m) {
				atoms[i] =atoms[i]+dr[m];
				double dp=angle(atoms[0],atoms[1],atoms[2]);
				atoms[i] = atoms[i]-2.0*dr[m];
				double dm=angle(atoms[0],atoms[1],atoms[2]);
				std::cout <<i<<"\t"<<m<<"\t"<<deriva[i][m]<<"\t"<<
						(cos(dp)-cos(dm))/0.0002<<std::endl;
				atoms[i] =atoms[i]+dr[m];
			}
		}
	std::cout << "dih 1-2-4:"<<
			torsion(atoms[0],atoms[1],atoms[2],atoms[3],&derivt) << std::endl;
	for(int i=0;i<4;++i){
			for(int m=0; m<3;++m) {
				atoms[i] =atoms[i]+dr[m];
				double dp=torsion(atoms[0],atoms[1],atoms[2],atoms[3]);
				atoms[i] = atoms[i]-2.0*dr[m];
				double dm=torsion(atoms[0],atoms[1],atoms[2],atoms[3]);
				std::cout <<i<<"\t"<<m<<"\t"<<derivt[i][m]<<"\t"<<
						(dp-dm)/0.0002<<std::endl;
				atoms[i] =atoms[i]+dr[m];
			}
		}
}


