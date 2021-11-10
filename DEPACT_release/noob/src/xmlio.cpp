/*
 * xmlio.cpp
 *
 *  Created on: 2018年8月9日
 *      Author: hyliu
 */
#include "xmlio.h"
#include <algorithm>
#include <cctype>
#include <locale>
using namespace subsitedesign;
 static inline std::string & trim(std::string &str){
	str.erase(str.begin(), std::find_if(str.begin(), str.end(), [](int ch) {
	        return !std::isspace(ch);}));
	str.erase(std::find_if(str.rbegin(), str.rend(), [](int ch) {
	        return !std::isspace(ch);
	    }).base(), str.end());
	return str;
}
static std::string cattextnodes(xmlpp::NodeSet nodes,bool dotrim=true){
	std::string str;
	for (auto n : nodes) {
		str = str + " "
				+ dynamic_cast<const xmlpp::TextNode*>(n)->get_content();
	}
	if(dotrim) trim(str);
	return str;
}
TargetStruct subsitedesign::xml2targetstruct(xmlpp::Node *ele){
	TargetStruct res;
	res.molname=cattextnodes(ele->find("molname/text()"));
	std::string atomtypestr=cattextnodes(ele->find("atomtypes/text()"));
	std::istringstream tsstr(atomtypestr);
	std::string at;
	while (tsstr.good()) {
		tsstr >>at;
		res.atomtypes.push_back(at);
	}
	std::string allrecords=cattextnodes(ele->find("pdbrecords/text()"));
	std::istringstream pdbsstr(allrecords);
	std::vector<NSPproteinrep::PdbRecord> records;
	std::string line;
	while(getline(pdbsstr,line)){
		if(line.substr(0,6) !="ATOM  " && line.substr(0,6) !="HETATM") continue;
		NSPproteinrep::PdbRecord nrecord(line);
		std::string e=std::string(nrecord.elementname);
		nrecord.atomname=trim(e)
				+std::to_string(nrecord.atomid);
		records.push_back(nrecord);
	}
	assert(records.size()==res.atomtypes.size());
	res.conformer=NSPproteinrep::make_aaconformer(records);
	return res;
}
TmpltSSAs::Alignment subsitedesign::xml2alignment(xmlpp::Node *ele) {
	TmpltSSAs::Alignment res;
	std::string pairs_str=cattextnodes(ele->find("alignedpairs/text()"));
	std::istringstream sstr(pairs_str);
	std::vector<int> iread;
	int ir;
	while (sstr.good()){
		sstr >> ir;
		iread.push_back(ir);
	}
	assert(iread.size() % 2 == 0);
	for (int i = 0; i < iread.size() / 2; ++i) {
		res.alignedpairs.insert(std::make_pair(iread[2 * i], iread[2 * i + 1]));
	}
	res.move3d=xml2move3d(ele->get_children("move3d").front());
	std::string rmsdstr=cattextnodes(ele->find("rmsd/text()"));
	sstr.str(rmsdstr);
	sstr.clear();
	sstr>>res.rmsd;
	return res;
}
void subsitedesign::addnodealignment(xmlpp::Node *parent, const TmpltSSAs::Alignment &algn) {
	xmlpp::Element *eleal=parent->add_child("alignment");
	xmlpp::Element *eleap=eleal->add_child("alignedpairs");
	std::ostringstream tsstr;
	for(auto &p:algn.alignedpairs){
		tsstr<<" "<<p.first<<" "<<p.second;
	}
	tsstr<<std::endl;
	eleap->set_child_text(tsstr.str());
	addnodemove3d(eleal,algn.move3d);
	xmlpp::Element *elermsd=eleal->add_child("rmsd");
	elermsd->set_child_text(std::to_string(algn.rmsd)+"\n");
}
Move3D subsitedesign::xml2move3d(xmlpp::Node *ele) {
	Move3D res;
	std::string trstr=cattextnodes(ele->find("translation/text()"));
	std::stringstream sstr(trstr);
	int idx = 0;
	while (sstr.good())
		sstr >> res.tvec[idx++];
	assert(idx==3);
	std::string rstr=cattextnodes(ele->find("rotation/text()"));
	sstr.str(rstr);
	sstr.clear();
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j)
			sstr >> res.rmatrix[i][j];
	}
	return res;
}
void subsitedesign::addnodemove3d(xmlpp::Node *parent,const Move3D& move3d){
	xmlpp::Element *elem3d=parent->add_child("move3d");
	xmlpp::Element *eletrans=elem3d->add_child("translation");
	xmlpp::Element *elert=elem3d->add_child("rotation");
	std::string trstr;
	for(auto d:move3d.tvec) trstr=trstr+" "+std::to_string(d);
	trstr=trstr+"\n";
	eletrans->add_child_text(trstr);
	std::ostringstream rstr;
	for(int i=0;i<3;++i){
		for(int j=0;j<3;++j) rstr<< " "<< move3d.rmatrix[i][j];
		rstr<<std::endl;
	}
	elert->add_child_text(rstr.str());
}
void subsitedesign::addnodetmpltssas(xmlpp::Node *parent, const TmpltSSAs & ssas){
	xmlpp::Element *ele=parent->add_child("tmpltssas");
	xmlpp::Element *elenmtp=ele->add_child("tmpltname");
	elenmtp->set_child_text(ssas.tmpltname);
	xmlpp::Element *elenmtg=ele->add_child("targetname");
	elenmtg->set_child_text(ssas.targetname);
	xmlpp::Element *elesp=ele->add_child("atomtypes");
	std::string spstr;
	for (auto & a:ssas.atomtypes){
		spstr=spstr+" "+a;
	}
	elesp->set_child_text(spstr);
	for(auto & a:ssas.alignments){
		addnodealignment(ele,a);
	}
}
TmpltSSAs subsitedesign::xml2tmpltssas(xmlpp::Node *ele){
	TmpltSSAs res;
	res.tmpltname=cattextnodes(ele->find("tmpltname/text()"));
	res.targetname=cattextnodes(ele->find("targetname/text()"));
	std::string spstr=cattextnodes(ele->find("atomtypes/text()"));
	std::istringstream sstr(spstr);
	int ir;
	std::string sr;
	while (sstr.good()){
		sstr >> sr;
		res.atomtypes.push_back(sr);
	}
	auto nodes=ele->find("alignment");
	for(auto &n:nodes){
		res.alignments.push_back(xml2alignment(n));
	}
	return res;
}
