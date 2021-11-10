/*
 * interaction.h
 *
 *  Created on: 2016年12月5日
 *      Author: hyliu
 */

#ifndef BACKBONE_LOOPINTERACTIONS_H_
#define BACKBONE_LOOPINTERACTIONS_H_
#include "backbone/backbonesite.h"

namespace NSPproteinrep {

class LoopInteractions {
public:
	typedef std::vector<BackBoneSite> Segment;
	struct Switches{
		bool steric_on;
		bool tetrasef_on;
		bool torsionmotif_on;
		bool torsion_on;
	};
	static Switches switches,switchstate_steric_torsion, switchstate_torsion_motif_tetra;
	static double steric_interaction(const Segment & seg);
	static double steric_interaction(const Segment & seg1, const Segment &seg2, bool successive=false);
	static double sef_interaction(const Segment & seg1, const Segment &seg2, bool successive=false);
	static double SCORECUT;
};

class SiteInteractions {
public:
	static double steric_interaction(const BackBoneSite & site);
	static double steric_interaction(const BackBoneSite & s1, const BackBoneSite & s2);
	static double steric_interaction(const BackBoneSite & s1, const BackBoneSite &s2, unsigned int sep);
};
double loop_context_steric_interaction(std::vector<std::vector<BackBoneSite>> *context,
		unsigned int headsegment, unsigned int tailsegment,
		const std::vector<BackBoneSite> & loop);
double loop_context_interaction(std::vector<std::vector<BackBoneSite>> *context,
		unsigned int headsegment, unsigned int tailsegment,
		const std::vector<BackBoneSite> & loop);
double loop_loop_sef(std::vector<BackBoneSite>::const_iterator & iter_a,
		int posi_a, int length_a,
		std::vector<BackBoneSite>::const_iterator & iter_b, int posi_b,int length_b);
}




#endif /* BACKBONE_LOOPINTERACTIONS_H_ */
