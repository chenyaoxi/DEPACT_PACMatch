/*
 * samplerefdist.cpp
 *
 *  Created on: 2018年8月23日
 *      Author: hyliu
 */

#include "analyzecontact.h"
#include "scorecontact.h"
#include "dstl/randomengine.h"
#include "geometry/calculators.h"
#include "geometry/rotation.h"

using namespace NSPgeometry;
using namespace myobcode;
using namespace subsitedesign;
void randommv(std::vector<XYZ> &obj,double rmax, double amax){
	auto &rng=NSPdstl::RandomEngine<>::getinstance().realrng();
	XYZ axis(rng,1.0);
	XYZ cent=center(obj);
	Rotation rot(QuaternionCrd(axis,rng()*amax),cent);
	XYZ trans(rng,rmax);
	for(auto &x:obj) {rot.apply(&x);x=x+trans;}
}
std::vector<XYZ> makeprobe(){
	std::vector<NSPgeometry::XYZ> probemol;
	probemol.push_back(XYZ( 0.0,0.0,0.0));
	probemol.push_back(XYZ(1.5,0.0,0.0));
	probemol.push_back(InternaltoXYZ(probemol[0],probemol[1],1.5,109.5*3.14159/180.0));
	probemol.push_back(InternaltoXYZ(probemol[0],probemol[1],
			probemol[2],1.5,109.5*3.14159/180.0,3.14159265));
	return probemol;
}
double clashenergy(const std::vector<std::string> &atypes, const std::vector<XYZ>&acrd,
		const std::vector<XYZ> &probe,double &rmin2,double &ene2){
	std::string NnOo("NnOo");
	double rcut2=5.0*5.0;
	double alpha=1.0;
	double ene1=0.0;
	ene2=0.0;
	rmin2=1000000.0;
	for(int i=0;i<atypes.size();++i){
		if(atypes.at(i)[0]=='H') continue;
		double radius=3.0;
		if(NnOo.find(atypes.at(i)[0])!=std::string::npos)radius=2.7;
		double rd2=radius*radius;
		XYZ xi=acrd[i];
		for(auto &p:probe){
			double rij2=(p-xi).squarednorm();
			if(rmin2>rij2) rmin2=rij2;
			double diff=(sqrt(rij2)-sqrt(rd2));
			if(diff<0){
				double eclash=exp((-diff)/0.3);
				double e1=eclash;
				if(eclash>4) e1=4;
				double e2=eclash;
				if(e2>15) e2=15;
				ene1+=e1;
				ene2+=e2;
			}
		}
	}
	ene2=ene2-ene1;
	if(rmin2>rcut2) ene1 += (rmin2-rcut2)*alpha;
	return ene1;
}

void collectdistr(const std::vector<std::string> &atypes, const std::vector<XYZ> &acrd,
		std::vector<XYZ> &probe,std::map<std::string,std::vector<double>> &distr1,
		std::map<std::string,std::vector<double>> &distr2,double ene2){
//	static DistBins bins({1.0,2.0,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.55,3.7,3.85,4.0,4.15,4.3});
	const DistBins &bins=DistBins::distbins();
	std::string NnOo("NnOo");
	double rmin2=1000000.0;
	std::map<std::string,std::vector<double>> tdist1;
	std::map<std::string,std::vector<double>> tdist2;
	tdist1["all"]=std::vector<double>(bins.num_bins(),0.0);
	tdist2["all"]=std::vector<double>(bins.num_bins(),0.0);
	for(int i=0;i<atypes.size();++i){
		if(atypes[i][0]=='H') continue;
		XYZ xi=acrd[i];
		int pidx=-1;
		double radius=3.0;
		if(NnOo.find(atypes.at(i)[0])!=std::string::npos)radius=2.7;
		for(auto &p:probe){
			++pidx;
			double rij2=(p-xi).squarednorm();
			if(rmin2>rij2) rmin2=rij2;
			double r=sqrt(rij2);
			if(r<bins.rmin()||r>bins.rmax()) continue;
			int bidx=bins.binid(r);
			std::map<std::string,std::vector<double>>*tdist=&tdist1;
			if(pidx==1||pidx==2) tdist=&tdist2;
			if(tdist->find(atypes[i])==tdist->end()){
				(*tdist)[atypes[i]]=std::vector<double>(bins.num_bins(),0.0);
			}
			double diff=r-radius;
			double ebias=0.0;
			double ediff=0.0;
			if(diff<0) {
				double e=exp(-diff/0.3);
				ebias=e;
				if(ebias>4) ebias=4;
				ediff=e;
				if(ediff>15) ediff=15;
			}
			ebias=-(ene2-ediff);
			double w=exp(ebias);
			(*tdist)[atypes[i]][bidx]+=w;
			(*tdist)["all"][bidx]+=w;
		}
	}
	if(rmin2<bins.rmin()*bins.rmin()||rmin2>bins.rmax()*bins.rmax()) return;
	double w=exp(rmin2-9.0);
	for(int i=0;i<2;++i){
		std::map<std::string,std::vector<double>>*tdist=&tdist1;
		std::map<std::string,std::vector<double>>*distr=&distr1;
		if(i==1) {tdist=&tdist2;distr=&distr2;}
		for(auto &td:*tdist){
			for(auto &d:td.second){
			d *=w;
			}
		}
		for(auto &td:*tdist){
			if(distr->find(td.first)==distr->end()) (*distr)[td.first]=std::vector<double>(bins.num_bins(),0.0);
			for(int i=0;i<bins.num_bins();++i){
				(*distr)[td.first][i] +=td.second[i];
			}
		}
	}
}
void myobcode::samplerefdist(const std::vector<std::string> &atypes, const std::vector<XYZ> &acrd,
		std::map<std::string,std::vector<double>> &distr1,
		std::map<std::string,std::vector<double>> &distr2){
	std::vector<XYZ> probe=makeprobe();
	auto &rng=NSPdstl::RandomEngine<>::getinstance().realrng(0.0,1.0);
	XYZ pcenter=center(probe);
	XYZ acenter=center(acrd);
	for(auto &x:probe) x=x+acenter-pcenter;
	int nsteps=50000;
	double rmin2;
	double e2;
	double eold=clashenergy(atypes,acrd,probe,rmin2,e2);
	eold+=rmin2;
	for(int n=0;n<nsteps;++n){
		std::vector<XYZ> probnew=probe;
		if(n>1000 &&(n%2==0)){
			collectdistr(atypes,acrd,probe,distr1,distr2,e2);
		}
		randommv(probnew,1.0,20.0);
		double e2new;
		double enew=clashenergy(atypes,acrd,probnew,rmin2,e2new);
		enew +=rmin2;
		bool accept=false;
		if(enew<eold) accept=true;
		else if(enew-eold<10.0) if(rng()<exp(eold-enew)) accept=true;
		if(accept) {
			probe=probnew;
			eold=enew;
			e2=e2new;
		}
	}
}



