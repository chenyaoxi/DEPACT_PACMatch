/*
 * backboneminimizer.h
 *
 *  Created on: 2017年8月11日
 *      Author: hyliu
 */

#ifndef BACKBONE_BACKBONEMINIMIZER_H_
#define BACKBONE_BACKBONEMINIMIZER_H_

#include "dstl/minimizer.h"
#include "backbone/backbonesite.h"
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
namespace NSPproteinrep{
class BackBoneMinimizer: public NSPdstl::SA_Minimizer<std::vector<BackBoneSite>> {
public:
	typedef std::vector<BackBoneSite> Chain;
	/**run energy minimization of a backbone chain
	 * @param control  parameters controlling the minimization protocol;
	 * @param chain  starting peptide chain
	 * @returns a pointer to the final chain after the run
	 *
	 */
	virtual std::shared_ptr<Chain> run(const NSPdstl::MinimizerControl & control, const Chain &chain);
	void writeminconf(const std::string & filename);
};
}
#endif /* BACKBONE_BACKBONEMINIMIZER_H_ */
