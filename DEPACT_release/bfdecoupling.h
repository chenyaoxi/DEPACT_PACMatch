/*
 * bfdecoupling.h
 *
 *  Created on: 2019年1月24日
 *      Author: yxchen
 */

#ifndef BFDECOUPLING_H_
#define BFDECOUPLING_H_

#include "basicfragment.h"
#include "twowaymatches.h"
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <memory>
using namespace subsitedesign;
namespace myobcode{
std::string bfdecoupling(
		const subsitedesign::BasicFragment *sf,OpenBabel::OBMol *mol);
}

#endif /* BFDECOUPLING_H_ */
