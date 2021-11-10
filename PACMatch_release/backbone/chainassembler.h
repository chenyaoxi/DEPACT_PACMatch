/*
 * chainassembler.h
 *
 *  Created on: 2017年1月4日
 *      Author: hyliu
 */

#ifndef BACKBONE_CHAINASSEMBLER_H_
#define BACKBONE_CHAINASSEMBLER_H_

#include "backbone/backbonesite.h"
#include <memory>
#include <map>
namespace NSPproteinrep {
/**
 * Add loops to link a set of fragments.
 * The fragments are store in SSelements_. The order of the fragments are stored in selements order_.
 * The generated loops are stored in loops_.
 * Note: because the loop closure algorithm changes the coordinates of the last residue of the
 * previous fragment and the first residue of the next fragment,
 * the first residue of a loop replaces the last residue of its N-terminal fragment,
 * and the last residue of a loop replaces the first residue of its C-terminal fragment.
 * In other words, each actual loop stored in loops_ will be two residue longer than the
 * length specified in looplengths.
 */
class ChainAssembler{
public:
	typedef std::vector<BackBoneSite> Fragment;
	struct Loop{
		std::shared_ptr<Fragment> loop{nullptr};
		unsigned int length{0};
		double score{10000000.0};
	};
	/**
	 * Addd fragment  one by one during setup.
	 */
	void addSSelement(const Fragment & ele);
	/**
	 * Access the elements
	 */
	const std::vector<Fragment> &SSelements() const {return SSelements_;}
	/**
	 * set order of SS element. Order is a list of integers. For example,set order to
	 * {1,0,2} coorespond to adding loops linking fragment 1 to 0, and 0 to 2, from the N
	 * to the C terminal.
	 */
	void setelementorder(const std::vector<unsigned int> & order);
	/**
	 * setup lengths of loop. The keys of the looplengths map are pairs of fragment IDs.
	 * For example, {{3,2},5} means a loop of length 5 linking the C terminal of fragment 3
	 * to N terminal of fragment 1.  Please keep in mind that the actual size of the backbone
	 * site vector storing the generated loop will be 5+2=7, with one extra residue at each
	 * end to replace the corresponding fragment residue.
	 */
	void setlooplengths(const std::map<std::pair<unsigned int,unsigned int>,unsigned int> & looplengths);

	/**
	 * rebuild the loop for the given head and tail fragment.
	 *
	 */
	bool rebuildloop(unsigned int headelement, unsigned int tailelement,unsigned int ntries);

	/**
	 * rebuild all loops connecting the fragments. This is usually called at the begining of
	 * a chain assembling procedure.
	 */
	bool rebuildallloops(unsigned int ntries);

	/**
	 * check if a loop has already been built.
	 */
	bool loopbuilt(unsigned int headelement, unsigned int tailelement) const {
		auto lp=loops_.find(std::make_pair(headelement,tailelement));
		if(lp==loops_.end()) return false;
		if(!lp->second.loop) return false;
		return true;
	}
	/**
	 * Returns a complete chain as a single fragment. It comprises input fragments and
	 * constructed loops in correct order. Tail and/or head residues of the fragments
	 * have been subsituted by cooreponding loop residues properly.
	 */
	Fragment getassembledchain(double *score);

private:
	std::vector<Fragment> SSelements_;
	std::map<std::pair<unsigned int, unsigned int>,Loop> loops_;
	std::vector<unsigned int>elementorder_;
	void makeloopcontext(unsigned int headelement,
			unsigned tailelement,std::vector<Fragment> *context,bool *complete);
};

/**
 * driver function to access Chain Assembler.
 * SS are the fragments or SS elements to be linker by loops.
 * order is the order of the elements.
 * loop lengths are the lengths of the loops (NOT include the boundary residues that
 * overlap the tail or head residue of the linked fragments)
 */
bool assembleSSelements(const std::vector<std::vector<BackBoneSite>> & SS,
		const std::vector<unsigned int> &order,
		const std::map<std::pair<unsigned int,unsigned int>,unsigned int> & looplengths);
}



#endif /* BACKBONE_CHAINASSEMBLER_H_ */
