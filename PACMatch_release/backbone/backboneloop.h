/*
 * backboneloop.h
 *
 *  Created on: 2016年12月7日
 *      Author: hyliu
 */

#ifndef BACKBONE_BACKBONELOOP_H_
#define BACKBONE_BACKBONELOOP_H_
#include "proteinrep/chaintree.h"
#include "backbone/backbonesite.h"
#include "pdbstatistics/phipsidistr.h"
#include "loopclosure/loopsolver.h"
#include "dstl/topn.h"

#include<memory>
namespace NSPproteinrep {
/**
 * The N and C terminal backbone coordinates for loop closure
 *
 */
struct FixEnds {
	NSPgeometry::XYZ r_c0; //needed for computing new phi1
	NSPgeometry::XYZ r_n1;
	NSPgeometry::XYZ r_ca1;
	NSPgeometry::XYZ r_c1;
	NSPgeometry::XYZ r_o1;
	NSPgeometry::XYZ r_ca3;
	NSPgeometry::XYZ r_c3;
	NSPgeometry::XYZ r_o3;
	NSPgeometry::XYZ r_n4; //needed for computing new psi3
	/**
	 * fixed coordinates needed by the tri-residue analytical loop closure solver
	 */
	std::vector<NSPgeometry::XYZ> fixedCrds() const {
		std::vector<NSPgeometry::XYZ> res;
		res.push_back(r_n1);
		res.push_back(r_ca1);
		res.push_back(r_ca3);
		res.push_back(r_c3);
		return res;
	}
};

/**
 * Deprecated.
 * One solution from the loop closure solver.
 * It re-wraps the solution from  loop closure solver into a single object.
 * Only the changed (phi,psi) pairs at the three rotated postions are stored in the created
 * object. Other (phi,psi) are contained in the input crd.
 */
struct ClosureSolution {
	unsigned int posi1, posi2, posi3;
	double phi1;
	double psi1;
	double phi2;
	double psi2;
	double phi3;
	double psi3;
	ClosureSolution(const FixEnds & ends, ChainTreeCrd *crd, unsigned int p1,
			unsigned int p2, unsigned int p3,
			const NSPloopclosure::LoopSolution &sol);
};
/**
 * It re-wraps the solution from  loop closure solver into a single object.
 * crdmap contains transformed atomic coordinates.
 * phi1 contains the new phi at the first position. This torsion angle is special because
 * it depends on, besides the loop atoms, the C atom from the residue right before the loop.
 */
struct CrdSolution {
	unsigned int posi1, posi2, posi3;
	double phi1;  //phi1 cannot be calculated from only the loop coordinates
	double psi3;  // It's more accurate to compute psi3 with the N crd of next residue
	std::map<typename ChainTreeTopology::AtomKey, NSPgeometry::XYZ> crdmap;
	CrdSolution(const FixEnds & ends, ChainTreeCrd *crd, unsigned int p1,
			unsigned int p2, unsigned int p3,
			const NSPloopclosure::LoopSolution &sol);
};
/**
 * This is the loop builder.
 */
class BackBoneLoop {
public:
	typedef typename ChainTreeTopology::AtomKey AtomKey;
//	BackBoneLoop() {
//		;
//	}
//	BackBoneLoop(const BackBoneSite &first, const BackBoneSite &last,
//			unsigned int length);
/**
 * context contains the fragments to be linked by the loop to be built.
 * context->at[headsegment] is right at the N-terminal of the loop.
 * context->at[tailsegment] is right at the C-terminal of the loop.
 *  The internal loop length will be two residue longer than the input
 *  length here, with the two extra boundary residues overlap
 *  the respective last(first) residue of the head(tail) segments.
 *  If loopsegment is given, it should be of the correct loop length.
 *  Then gly and pro will not be chosen randomly. Their locations in input segments will
 *  be retained when build the loop.
 */
	BackBoneLoop(std::vector<std::vector<BackBoneSite>> *context,
			unsigned int headsegment, unsigned int tailsegment,
			unsigned int length, const std::vector<BackBoneSite> & loopsegment =
					std::vector<BackBoneSite>());
/**
 * Initialize internal coordinates of the loop.
 * (phi,psi) of the first and last residues are taken from first and last.
 * Coordinates of the N-terminal loop atoms are taken from the input first residue.
 */
	void initConf(const BackBoneSite &first, const BackBoneSite &last);
/**
 * assign gly and pro type to positions
 * These residue types lead to speical torsional preferences.
 * Pro may also allow trans peptide bond.
 */
	void chooseglypro();

/**
 * randomize the backbone torsional angles.
 * (phi,psi) are taken from respective distributions statistically generated from database.
 */
	void randomConf(bool choosegp = false);

/**
 * randomize a portion of the torsional angles
 */
	void partial_randomConf();
/**
 * scores a single loop closure solution.
 * Only the changing (phi,psi) pairs of the three residues are scored in this method.
 */
	double scoreSolution(const ClosureSolution & s) const;
/**
 * scores a single loop closure solution.
 * Besides torsion score, it also calls loop_context_interaction to score the interaction of the built loop with
 * context fragments.
 */
	double scoreCrdSolution(const CrdSolution &s);
/**
 * calculate torsion score. Each pair of (phi,psi) are scored independently
 * based on the respective statistical distribution,and then summed.
 */
	double torsionscore() const;
/**
 * Obtain loop closure solutions
 * posi is the location of middle residue in the three-residue
 * loop closure algorithm.
 */
	void closureSolutions(unsigned int posi,
			std::vector<CrdSolution> &solns) const;
/**
 * Obtain loop closure solutions
 */
	void closureSolutions(unsigned int p1, unsigned int p2, unsigned int p3,
			std::vector<ClosureSolution> &solns) const;

/**
 * Choose from a set of solutions based on a Boltzmann distribution calculated from
 * the scores.
 */
	double chooseApplySolution(std::vector<CrdSolution> &solns,
			double temperature);
/**
 * generate a new loop conformation.
 * It first randomize the loop conformation,then obtain loop closure solutions,
 * The call chooseApplySolution to select a single solution.
 * Returns false if failed.
 */
	bool newLoopConf(std::vector<BackBoneSite> *loop, double *sc);
/**
 * calls newLoopConf repetitively (at least ncopies*ncandidates times)
 * to generate ncopies of loop conformations. The best ncopies of built loops are
 * copied to the container position started from the input iterator loopiter.
 * If scores pointer is the corresponding loop scores are delivered
 */
	unsigned int getLoops(unsigned int ncopies, unsigned int ncandidates,
			typename std::vector<std::vector<BackBoneSite>>::iterator loopiter,std::vector<double> *scores=nullptr);
/**
 * Sample a low energy loop. Try to generate as many as ncandidates*nsamples allowed loop conformations (successfully
 * closed, clash-free), from which the top ncandidates with the lowest torsion scores will be selected.
 * Each of the candidates is then evaluated by the total score (including torsion, motif and tetrasef scores).
 * The one with the lowest total energy is delivered to resultloop, its total score returned.
 */
	double sampleALoop(unsigned int ncandidates,unsigned int nsamples,
			std::unique_ptr<std::vector<BackBoneSite>> resultloop);
/**
 * A helper to obtain the phi-psi distribution statistics
 */
	const NSPpdbstatistics::PhiPsiDistr *phipsidistr(unsigned int posi) const;
/**
 * Turn built loop into backbonesite objects
 */
	BackBoneSite makeSite(unsigned int posi, unsigned int shift = 1,
			std::string pdbid = "LLP", char chainid = 'A', std::string resname =
					"ALA") const;
	void makeAllsites(std::vector<BackBoneSite> *loop,
			unsigned int shift) const {
		for (unsigned int posi = 0; posi < length_; ++posi) {
			loop->push_back(makeSite(posi, shift));
		}
	}
	void setsaveacceptable(double refscore){
		refscore_=refscore;
		saveacceptable_=true;
	}
	int saveTopSolutions(std::vector<CrdSolution> &solns);
	double copyNScoreLoop(std::vector<BackBoneSite>::const_iterator &iter);
	double loopScore();
	int batchFindSolutions(unsigned int ncopies);
	void clearTopSolutions() {topsolutions_.clear();}
private:
	std::shared_ptr<ChainTreeTopology> chaintopo_;
	unsigned int length_;
	std::shared_ptr<ChainTreeCrd> crd_;
	std::vector<bool> isgly_;
	std::vector<bool> ispro_;
	FixEnds ends_;
	BackBoneSite first_, last_;
	std::vector<std::vector<BackBoneSite>> *contextsegments_ { nullptr };
	unsigned int headsegment_;
	unsigned int tailsegment_;
	bool chooseglypro_{true};
	bool firstcispro_{false};
	bool saveacceptable_{false};
	double refscore_{1.e20};
	unsigned int nsolns_save{100};
	NSPdstl::TopN<std::shared_ptr<std::vector<BackBoneSite>>> topsolutions_;
	std::vector<std::shared_ptr<std::vector<BackBoneSite>>>acceptablesolutions_;
	std::vector<double> acceptablescores_;
};

}

#endif /* BACKBONE_BACKBONELOOP_H_ */
