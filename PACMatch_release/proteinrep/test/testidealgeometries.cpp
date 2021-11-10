/*
 * testidealgeometries.cpp
 *
 *  Created on: 2016年11月17日
 *      Author: hyliu
 */

#include "proteinrep/idealgeometries.h"
#include <vector>
#include <iostream>
using namespace NSPproteinrep;
int main() {
	const std::string filename{"idealgeometries.dat"};
	IdealGeometries & igdat=IdealGeometries::getGlobalInstance(filename);
	double deg=igdat.degree;
	std::vector<typename IdealGeometries::Bond> bonds{{"N","H"},{"H","N"},{"N","CA"},
		{"C","O"},{"N","CA"},{"CA","C"},{"CH3","N"},{"CH3","C"}
	};
	std::vector<typename IdealGeometries::Angle> angles{{"C","N","H"},{"CA","C","O"},{"N","CA","C"},
	};
	std::vector<typename IdealGeometries::RelativeTorsion> rtorsions{{"CA","C","O","N"},
		{"C","N","CA","H"},{"N","CA","CB","C"}};

	for(auto &b: bonds) {
		std::cout <<b.first <<"-" <<b.second << ": " <<igdat.idealLength(b.first,b.second)<<std::endl;
	}
	for(auto &a: angles) {
		std::cout <<a.ia <<"-" <<a.ja <<"-"<<a.ka<< ": " <<igdat.idealAngle(a.ia,a.ja,a.ka)/deg<<std::endl;
	}
	for(auto &rt:rtorsions){
		std::cout <<rt.ja<<"-" <<rt.ka<<"-"<<rt.la1 <<"&"<<rt.la2<<": "
				<<igdat.idealRelativeTorsion(rt.ja,rt.ka,rt.la1,rt.la2)/deg <<std::endl;
	}
}

