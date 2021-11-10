/*
 * clusterpepconformer.cpp
 *
 *  Created on: 2017年6月29日
 *      Author: hyliu
 */

#include "backbone/pepconformer.h"
#include "kmeans/runkmeans.h"
using namespace NSPproteinrep;
using namespace NSPgeometry;
using namespace NSPkmeans;

int main(int argc,char **argv){
	std::vector<BackBoneSite> sites;
	readbackbonesites(std::string(argv[1]), sites);
	std::vector<std::vector<double>> data;
	std::vector<int> positions;
//	std::vector<PepConformer> conformers;
	int maxconformers=10000;
	int peplength=5;
	double deg=3.14159265/180.0;
	for(int posi=0; posi<=sites.size()-peplength;++posi){
		if(!fragstartsite(sites.begin()+posi,sites.end(), peplength,std::string(),false))
			continue;
//		conformers.push_back(make_conformer(peplength,sites,posi+peplength/2));
		positions.push_back(posi);
		data.push_back(std::vector<double>());
		std::vector<double> & vec=data.back();
		double atemp;
		atemp=sites[posi].psi()*deg;
		vec.push_back(cos(atemp));
		vec.push_back(sin(atemp));
		atemp=sites[posi].omiga()*deg;
		vec.push_back(cos(atemp));
		for(int i=1;i<peplength-1; ++i) {
			atemp=sites[posi+i].phi()*deg;
			vec.push_back(cos(atemp));
			vec.push_back(sin(atemp));
			atemp=sites[posi+i].psi()*deg;
			vec.push_back(cos(atemp));
			vec.push_back(sin(atemp));
			atemp=sites[posi+i].omiga()*deg;
			vec.push_back(cos(atemp));
		}
		atemp=sites[posi+peplength-1].phi()*deg;
		vec.push_back(cos(atemp));
		vec.push_back(sin(atemp));
//		if(conformers.size() >= maxconformers) break;
		if(positions.size()>= maxconformers) break;
	}
/*	std::vector<double> ra(4*peplength,0.0);
	for(auto &c:conformers){
		std::vector<XYZ> & lcrd=c.localcrd();
		for(int i=0; i<4*peplength;++i) {
			ra[i] += sqrt(lcrd[i].squarednorm());
		}
	}
	for(int i=0; i<4*peplength;++i) {
		ra[i] = ra[i]/(double) conformers.size();
		if(ra[i]<0.01) ra[i]=1.0;
	}*/
	NSPkmeans::RunKMeans rkmeans;
	int dim=data[0].size();
	int maxpoints=positions.size();
	rkmeans.init_data(dim,maxpoints);
/*
	for(auto &c:conformers){
		std::vector<double> point;
		std::vector<XYZ> & lcrd=c.localcrd();
		int idx=0;
		for(auto &atm:lcrd) {
			double w=1.0/ra[idx++];
			point.push_back(atm.x_*w);
			point.push_back(atm.y_*w);
			point.push_back(atm.z_*w);
		}
		rkmeans.addpoint(point);
	}*/
	for(auto &d:data){
		rkmeans.addpoint(d);
	}
	rkmeans.finishdata();
//	for(int k=5; k<40; ++k) {
	int k=25;
		rkmeans.setrun(k,NSPkmeans::RunKMeans::HYBRID);
		rkmeans.run();
		std::cout<<"cluster " << k<<'\t'<<rkmeans.getdistortion()<<std::endl;
		auto &clusters=rkmeans.clusters();
		int clusterid=0;
		for(auto &c:clusters){
				std::cout<<c.distortion;
				std::cout<<'\t'<<c.members.size()<<std::endl;
				for(int i=0; i<dim;++i)
				std::cout <<'\t'<< c.center[i];
				std::cout<<std::endl;
				std::ofstream ofs;
				std::string filename="cluster_" + std::to_string(clusterid++) +".torsions";
				ofs.open(filename.c_str());
				int mid=0;
				for(auto & idx_d:c.members){
					int idx=idx_d.first;
/*					std::vector<BackBoneSite> csites;
					conformers[idx].tobackbonesites(csites,conformers[idx].localcrd());
					ofs<<"MODEL  " <<++mid<<std::endl;
					writeSitesToPDB(ofs,csites);
					ofs<< "ENDMDL" <<std::endl;*/
					int p=positions[idx];
					ofs <<sites[p].psi()<<' '<< sites[p].omiga();
					for(int i=1;i<peplength-1;++i) {
						ofs<<' '<<sites[p+i].phi()<<' '<<sites[p+i].psi()<<' '<<sites[p+i].omiga();
					}
					ofs <<' '<<sites[p+peplength-1].phi()<<std::endl;
				}
				ofs.close();
			}
//	}
}
