/*
 * myobxmlio.cpp
 *
 *  Created on: 2018年8月12日
 *      Author: hyliu
 */

#include "myobxmlio.h"

using namespace xmlpp;
using namespace OpenBabel;
using namespace myobcode;
void myobcode::addnodetargetstruct(xmlpp::Node *parent,
	OpenBabel::OBMol *mol,const std::vector<std::string> &atomtypes){
	Element *elestr=parent->add_child("targetstruct");
	Element *eletitle=elestr->add_child("molname");
	eletitle->set_child_text(std::string(mol->GetTitle()));
	Element *eleatomtypes=elestr->add_child("atomtypes");
	Element *elepdbrecords=elestr->add_child("pdbrecords");
	std::ostringstream tstr;
	for(auto &at:atomtypes){
		tstr<<" "<<at;
	}
	eleatomtypes->set_child_text(tstr.str());
	std::ostringstream pdbstr;
	OBConversion conv;
	conv.SetOutFormat("pdb");
	conv.Write(mol,&pdbstr);
	elepdbrecords->set_child_text(pdbstr.str());
}


