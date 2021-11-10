/*
 * scorepep.h
 *
 *  Created on: 2017年6月30日
 *      Author: hyliu
 */

#ifndef BACKBONE_SCOREPEP_H_
#define BACKBONE_SCOREPEP_H_
#include "backbone/backbonesite.h"
#include "backbone/torsionvectortree.h"
#include <set>
namespace NSPproteinrep {
class ScorePep {
public:
	enum conftype{COIL,HELIX,STRAND,ALL};
	void buildtree(const std::vector<BackBoneSite> &sites,int peplength);
	void buildreftree(long npoints);
	void buildrefphipsitrees(long npoints);
	double score(const std::vector<double> & conf,
			double *srx,double *sr_ana,double *d2m,double *d20);
	double neighborsum(const std::vector<double> & conf,
			const std::vector<double>&weights,TorsionVectorTree &tree);
	double refneighborsum(const std::vector<double> &conf);
	std::vector<double> extractconf(const std::vector<BackBoneSite> &sites, long posi,int length);
	std::vector<double> sampleconf(int l,int peptype=ALL);
	long templatesize() const {return templateconfs_->size();}
	long refsize() const {return refconfs_->size();}
	TorsionVectorTree &tree() {return tree_;}
	TorsionVectorTree &reftree(){ return reftree_;}
	double refprobability(const std::vector<double> &conf) const;
	void addexcludedpb_(const std::vector<char> &exclude={'d','m'}){
		for(auto pb:exclude) excludedpb_.insert(pb);
	}
	ScorePep(int conftype=COIL):conftype_(conftype){
		if(conftype_==COIL) addexcludedpb_();
	}
private:
	int length_{5};
	int conftype_{COIL};
	double pbmweight_{0.2};
	std::set<char> excludedpb_;
	TorsionVectorTree tree_;
	TorsionVectorTree reftree_;
	TorsionVectorTree refphi_;
	TorsionVectorTree refpsi_;
	TorsionVectorTree refphipsi_;
	TorsionVectorTree refphipsip_;
	TorsionVectorTree refphipsin_;
	std::shared_ptr<std::vector<std::vector<double>>> templateconfs_;
	std::shared_ptr<std::vector<double>> templateweights_;
	double templatewtot_{0.0};
	std::shared_ptr<std::vector<std::vector<double>>> refconfs_;
	std::shared_ptr<std::vector<std::vector<double>>> refphis_;
	std::shared_ptr<std::vector<std::vector<double>>> refpsis_;
	std::shared_ptr<std::vector<std::vector<double>>> refphipsis_;
	std::shared_ptr<std::vector<std::vector<double>>> refphipsips_;
	std::shared_ptr<std::vector<std::vector<double>>> refphipsins_;
	std::shared_ptr<std::vector<double>> refweights_;
	double refwtot_{0.0};
	static double pairscore(const std::vector<double> &conf, const std::vector<double> & tmpl);
};

}




#endif /* BACKBONE_SCOREPEP_H_ */
