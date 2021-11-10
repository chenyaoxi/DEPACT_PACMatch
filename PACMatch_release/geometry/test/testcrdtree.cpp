/*
 * testcrdtree.cpp
 *
 *  Created on: 2016年11月15日
 *      Author: hyliu
 */

#include "geometry/crdtree.h"
#include "geometry/rotation.h"
#include <iostream>
using namespace NSPgeometry;
int main() {
	typedef TopoTree<std::string> Topology;
	typedef CrdTree<std::string> Coord;
	typedef typename TopoTree<std::string>::Tree Tree;
	typedef typename Tree::Pointer TreePointer;

	std::vector<std::string> atomkeys{"ME1","C1","O1","N1","H1","CA","CB","C2","O2","N2","H2","ME2"};
	Topology topo(atomkeys[0]); //ME1
	Coord crd(&topo);
	double phi=-60.0;
	double psi=120.0;
	double deg=3.14159265/180.0;
	crd.addAtom(atomkeys[0],IntCrd());
	topo.attachAtom(atomkeys[0],atomkeys[1],true); //c1
	crd.addAtom(atomkeys[1],IntCrd(1.53));
	topo.attachAtom(atomkeys[1],atomkeys[2]); //O1
	crd.addAtom(atomkeys[2],IntCrd(1.23,121.0*deg));
	topo.attachAtom(atomkeys[1],atomkeys[3],true); //N1
	crd.addAtom(atomkeys[3],IntCrd(1.33,115.0*deg));
	topo.attachAtom(atomkeys[3],atomkeys[4]); //H1
	crd.addAtom(atomkeys[4],IntCrd(1.0,123.0*deg,0.0));
	topo.attachAtom(atomkeys[3],atomkeys[5],true); //CA
	crd.addAtom(atomkeys[5],IntCrd(1.47,122.0*deg,180.0*deg));
	topo.attachAtom(atomkeys[5],atomkeys[6]); //CB
	crd.addAtom(atomkeys[6],IntCrd(1.53,109.5*deg,(phi+130)*deg));
	topo.attachAtom(atomkeys[5],atomkeys[7],true); //C2
	crd.addAtom(atomkeys[7],IntCrd(1.53,109.5*deg, phi*deg));
	topo.attachAtom(atomkeys[7],atomkeys[8]); //O2
	crd.addAtom(atomkeys[8],IntCrd(1.23,121.0*deg, (psi+180.0)*deg));
	topo.attachAtom(atomkeys[7],atomkeys[9],true); //N2
	crd.addAtom(atomkeys[9],IntCrd(1.33,115.0*deg, psi*deg));
	topo.attachAtom(atomkeys[9],atomkeys[10]); //H2
	crd.addAtom(atomkeys[10],IntCrd(1.0,123.0*deg, 0.0));
	topo.attachAtom(atomkeys[9],atomkeys[11],true); //ME2
	crd.addAtom(atomkeys[11],IntCrd(1.47,122.0*deg,180.0*deg));

	std::map<std::string, XYZ> crd0;
	crd0.insert(std::make_pair(atomkeys[0],XYZ()));  //ME1
	crd0.insert(std::make_pair(atomkeys[1],InternaltoXYZ(crd0.at(atomkeys[0]),1.53))); //C1
	crd0.insert(std::make_pair(atomkeys[2],
			InternaltoXYZ(crd0.at(atomkeys[1]),crd0.at(atomkeys[0]),1.23,121.0*deg))); //O1
	XYZ temp=InternaltoXYZ(crd0.at(atomkeys[1]),crd0.at(atomkeys[0]),1.33,115.0*deg);
	Rotation r(QuaternionCrd(crd0[atomkeys[1]]-crd0[atomkeys[0]],180.0),crd0[atomkeys[1]]);
	crd0.insert(std::make_pair(atomkeys[3],r.applytoCopy(temp))); //N1

	crd.copyCrdMap(crd0);
//	crd.calcXYZ();

	std::map<std::string,XYZ> &crdmap=crd.crdmap(true);
	for(auto k: atomkeys) {
		std::cout << k<<"\t"<<crdmap[k].toString() <<std::endl;
	}

	Coord crd1(&topo);
/*
	for(auto & c:crd0) {
		c.second = c.second +XYZ(5,5,5);
	}
	*/
	crd1.copyCrdMap(crd0);
	crd1.updateInternal(crdmap);
	crd1.rotateBond(atomkeys[5],-60*deg);
	crd1.rotateBond(atomkeys[7],30*deg);
	crdmap=crd1.crdmap(true);
	for(auto k: atomkeys) {
		std::cout << k<<"\t"<<crdmap[k].toString() <<std::endl;
	}
	crd.copyCrdMap(crdmap);
	crd.calcInternal();
	std::map<std::string,IntCrd> &icrdmap=crd.intcrdmap();
	for(auto k: atomkeys) {
			typename Topology::AIJK &aijk=topo.aijk(k);
			IntCrd &ic=icrdmap[k];
			std::cout << k<<"\t" << aijk.ka <<"\t" <<ic.distance <<"\t"<<aijk.ja<<"\t" << ic.angle/deg
					<<"\t"<<aijk.ia <<"\t"<<ic.torsion/deg
					<<std::endl;
		}

}


