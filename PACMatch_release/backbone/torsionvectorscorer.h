/*
 * torsionvectorscorer.h
 *
 *  Created on: 2016年12月28日
 *      Author: hyliu
 */

#ifndef BACKBONE_TORSIONVECTORSCORER_H_
#define BACKBONE_TORSIONVECTORSCORER_H_
#include "backbone/torsionvectortree.h"
#include "backbone/backbonesite.h"

namespace NSPproteinrep {

class TorsionVectorScorer {
public:
	const static double motifunfoundscore;
	const static double d2max;
	const static double d20;
	const static double d2sigma;
	const static double resolution;
	static TorsionVectorScorer &getinstance(std::vector<BackBoneSite> *sites=nullptr) {
		static TorsionVectorScorer instance;
		static const int LENGTH{4};
		static bool initialized{false};
		if(!initialized){
			assert(sites != nullptr);
			instance.init(sites,LENGTH);
			initialized=true;
		}
		return instance;
	}
	void init(std::vector<BackBoneSite> * sites, int length);
	double score(const std::string &motifname,const std::vector<double> & torsions);
	double score(std::vector<BackBoneSite>::iterator begin);
	double scorerange(std::vector<BackBoneSite>::iterator begin,std::vector<BackBoneSite>::iterator end);
	int length()const {return length_;}
private:
	std::map<std::string, TorsionVectorTree> templatetrees_;
	std::shared_ptr<std::vector<std::vector<double>>> templatevectors_;
	int length_;
	double neighborcut2(){
		return d2max*2.0*(double) length_;
	}
	double score(const std::vector<double> & torsions,const std::vector<double> &tmpl);
};

template<typename ITER>
std::string getmotifname(ITER iter,int length){
		std::string motif;
		std::string omigaseq;
		for(int i=0;i<length; ++i){
			if(iter->resname=="GLY") motif=motif+"G";
			else if(iter->resname=="PRO")motif=motif+"P";
			else motif=motif+"X";
			if(iter->omiga() <-90|| iter->omiga()>90){
				omigaseq +="c";
			} else {
				omigaseq +="t";
			}
			++iter;
		}
		return motif+omigaseq;
}

}



#endif /* BACKBONE_TORSIONVECTORSCORER_H_ */
