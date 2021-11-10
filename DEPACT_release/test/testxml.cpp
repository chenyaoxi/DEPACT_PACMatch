/*
 * testxml.cpp
 *
 *  Created on: 2018年8月7日
 *      Author: hyliu
 */

#include <libxml++/libxml++.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
using namespace OpenBabel;
int main(int argc,char **argv){
	std::string initdoc("<molwrapper></molwrapper>");
	OBConversion obconversion;
	OBMol mol;
	obconversion.SetInAndOutFormats("sdf","sdf");
	std::string targetsdf(argv[1]);
	std::cout <<argv[1]<<' '<<targetsdf<<std::endl;
	std::ifstream ifs(argv[1]);
	obconversion.Read(&mol,&ifs);
	std::stringstream oss;
	obconversion.Write(&mol,&oss);
//	bool notatend = obconversion.ReadFile(&mol,"target.sdf");
	std::string moltitle(mol.GetTitle());
	xmlpp::DomParser parser;
	parser.parse_memory(initdoc);
	xmlpp::Node * root=parser.get_document()->get_root_node();
	xmlpp::Element *molnode=root->add_child("MolDATA");
	molnode->set_attribute("Title",moltitle);
	molnode->set_attribute("format","sdf");
	molnode->add_child_cdata(oss.str());
	xmlpp::Element *kw=dynamic_cast<xmlpp::Element *>(root->get_first_child("keywords"));
/*	for(auto & w: keywords){

	}*/
	auto nodes=root->find("//keywords/text()[2]");
	std::cout<<kw->get_path()<<std::endl;
	for(xmlpp::Node* n:nodes){
		const auto data=dynamic_cast<const xmlpp::CdataNode*>(n);
		if(data) std::cout<< data->get_content()<<std::endl;
	}
	parser.get_document()->write_to_stream_formatted(std::cout);
}


