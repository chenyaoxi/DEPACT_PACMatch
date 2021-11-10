/*
 * myobxmlio.h
 *
 *  Created on: 2018年8月12日
 *      Author: hyliu
 */

#ifndef MYOBXMLIO_H_
#define MYOBXMLIO_H_
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <libxml++/libxml++.h>
#include <sstream>
#include <cassert>
namespace myobcode{
void addnodetargetstruct(xmlpp::Node *parent, OpenBabel::OBMol *mol,
		const std::vector<std::string> &atomtypes);
}



#endif /* MYOBXMLIO_H_ */
