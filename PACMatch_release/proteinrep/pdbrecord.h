/*
 * pdbatomrecord.h
 *
 *  Created on: 2016年11月16日
 *      Author: hyliu
 */

#ifndef PROTEINREP_PDBRECORD_H_
#define PROTEINREP_PDBRECORD_H_
#include <string>
#include <vector>
namespace NSPproteinrep{
#ifndef NOID
#define NOID -100000
#endif
#ifndef CHAINID_DEFAULT
#define CHAINID_DEFAULT 'A'
#endif
/*! data in an ATOM or HETATM record
 *
 */
struct PdbRecord {
	/*!fields in one line
	 *
	 */
	enum Field {
		LABEL,
		ATOMID,
		NAMESYMBOL,
		NAMEMODIFIER,
		CONFORMERID,
		RESIDUENAME,
		CHAINID,
		RESIDUEID,
		INSERTIONID,
		X,
		Y,
		Z,
		OCCUPATION,
		BFACTOR,
		SEGMENT,
		ELEMENTNAME
	};
	/*!default constructor
	 *
	 */
	PdbRecord() {
		;
	}
	PdbRecord(const std::string & line);
	void init(const std::vector<std::string> &fields);
	std::string toString() const;
	std::string label { "" };
	std::string atomname { "" };  //atomname=namesymbol+namemodifier
	std::string namesymbol { "" };
	std::string namemodifier { "" };
	std::string residuename { "" };
	int atomid { NOID };
	char chainid { ' ' };
	int residueid { NOID };
	char conformerid { ' ' };
	char insertionid { ' ' };
	double x { 0.0 }, y { 0.0 }, z { 0.0 };
	double occupation { 1.0 };
	double bfactor { 0.0 };
	char elementname[2] { ' ', 'X' };
	bool failure {false};
};
template<typename ATOMKEY,typename XYZ>
PdbRecord make_pdbrecord(const typename ATOMKEY::Key & key, const XYZ & xyz,int atomid=0, int residueidshift=0,
		int elementsymbolsize=1) {
	PdbRecord record;
	const std::vector<char> chainids{'A','B','C','D','E','F','G','H','I','J'};
	record.label="ATOM";
	record.atomname=ATOMKEY::atomName(key);
	record.namesymbol=record.atomname.substr(0,elementsymbolsize);
	if(elementsymbolsize > 1) {
		record.elementname[0]= record.namesymbol[0];
		record.elementname[1]=record.namesymbol[1];
	} else {
		record.elementname[1]=record.namesymbol[0];
	}
	record.namemodifier=record.atomname.substr(elementsymbolsize);
	record.residuename=ATOMKEY::residueName(key);
	record.atomid=atomid;
	record.residueid = ATOMKEY::posiNumber(key) + residueidshift;
	record.x=xyz[0];
	record.y=xyz[1];
	record.z=xyz[2];
	int chainid=ATOMKEY::chainNumber(key);
	if(chainid <10) record.chainid=chainids[chainid];
	return PdbRecord(record.toString());
}
}

#endif /* PROTEINREP_PDBRECORD_H_ */
