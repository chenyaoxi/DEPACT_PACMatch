/*
 * closingloop.h
 *
 *  Created on: 2017年4月23日
 *      Author: hyliu
 */

#ifndef BACKBONE_CLOSINGLOOP_H_
#define BACKBONE_CLOSINGLOOP_H_
#include "backbone/backbonesite.h"
#include "loopclosure/loopsolver.h"
#include <memory>

namespace NSPproteinrep {
class MainChain;
class ClosingLoop {
public:
	enum RegenMode {NEW,MUTATE,KEEP};
	enum ModeType {NEWSEQ_NEWCONF,NEWSEQ_MUTATECONF,MUTATESEQ_CONF,KEEPSEQ_NEWCONF,
		KEEPSEQ_MUTATECONF};
	typedef std::vector<BackBoneSite> Solution;
	ClosingLoop(){;}
	ClosingLoop (const MainChain & host,int startposi,int endposi,
			int modestype,int pertposi=-1);
	int solve();
	bool unsolvable() {return allowedp2_.empty();}
	void clearsolutions(){solutions_.clear();p2ofsolutions_.clear();}
	void getsitessolution(int i,std::vector<BackBoneSite> *result);
private:
	MainChain tmplt_;
	bool nextpro_{false};  //needed for choosing phipsidistr of last site.
	double lastomiga_{180.0};
	std::vector<NSPloopclosure::LoopSolution> solutions_;
	std::vector<int> p2ofsolutions_;
	std::vector<int> allowedp2_; //middle positions of non-rigid local conformation
	std::vector<std::string> seq_;
	std::vector<BackBoneSite> refsites_;
	std::vector<NSPgeometry::XYZ> fixcrds_;
	NSPgeometry::XYZ rc0_; //needed to calculate phi of first loop site
	NSPgeometry::XYZ rn4_; //needed to calculate psi of last loop site;

	RegenMode seqmode_{NEW};
	RegenMode confmode_{NEW};
	int pertposi_{-1};
	int selectedpertposi_{-1};
	void setmodes(int modestype);
	void initseq();
	void initconf();
	void initseqconf();
	int addsolutions(int posi);
};
std::string chooseresiduetype();
void choosephipsi(std::string restype,std::string nextresiduetype,double *phi,double *psi,bool *cispro);
}

#endif /* BACKBONE_CLOSINGLOOP_H_ */
