/*
 * chainpackminimizer.h
 *
 *  Created on: 2017年8月17日
 *      Author: hyliu
 */

#ifndef BACKBONE_CHAINPACKMINIMIZER_H_
#define BACKBONE_CHAINPACKMINIMIZER_H_


#include "dstl/minimizer.h"
#include "backbone/chainpack.h"
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
namespace NSPproteinrep{
class ChainPackMinimizer: public NSPdstl::SA_Minimizer<std::vector<std::vector<BackBoneSite>>> {
public:
	typedef std::vector<BackBoneSite> Chain;
	virtual std::shared_ptr<std::vector<std::vector<BackBoneSite>>> run(const NSPdstl::MinimizerControl & control,
			const std::vector<std::vector<BackBoneSite>> &chains);
	void writeminconf(const std::string & filename);
};
}


#endif /* BACKBONE_CHAINPACKMINIMIZER_H_ */
