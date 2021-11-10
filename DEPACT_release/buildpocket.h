/*
 * buildpocket.h
 *
 *  Created on: 2019年7月3日
 *      Author: yxchen
 */

#ifndef BUILDPOCKET_H_
#define BUILDPOCKET_H_

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>
#include "basicfragment.h"
#include "mmmatches.h"
#include "atomtypessmarts.h"
#include "tmpltssas.h"
#include "cliquer/cliquer.h"
#include "subsite.h"
#include "proteinrep/aaconformer.h"
#include <openbabel/fingerprint.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
using namespace myobcode;
using namespace subsitedesign;
using namespace OpenBabel;
using namespace NSPproteinrep;
using namespace std;
#define DIST2CLASH 25.0

namespace buildpocket
{
std::map<std::string, std::string> readpar(std::string parname);
std::string getp(std::map<std::string, std::string> parmap, std::string name);
bool similarity(OBMol targetmol, OBMol templatemol, double perc);
double similar(OBMol targetmol, OBMol templatemol);
std::vector<OBMol> buildtmpmols(std::vector<int> bfneeds,
		std::vector<std::vector<std::string>> bf_names,
		OBConversion obconversion, OBMol targetmol, double perc);
std::vector<std::vector<OBMol>> buildtmpmols_sep(std::vector<int> bfneeds,
		std::vector<std::vector<std::string>> bf_names,
		OBConversion obconversion, OBMol targetmol, double perc);
void makessites(TmpltSSAs ssa, TargetStruct tgt, double weight,
		std::vector<std::shared_ptr<Subsite>> &test_ssites, int &idx,
		map<int, vector<string>> &covernum_codes,
		map<string, map<double, vector<int>>> &code_ids);
//void unredundent(std::vector<std::shared_ptr<Subsite>> pre_ssites,
//		std::vector<std::shared_ptr<Subsite>> &ussites,
//		map<string, map<double, vector<int>>> pcode_ids,
//		map<string, map<double, vector<int>>> &ucode_ids, double rmsd);
void topssites(std::vector<std::shared_ptr<Subsite>> pre_ssites,
		std::vector<std::shared_ptr<Subsite>> &ssites,
		map<string, map<double, vector<int>>> pcode_ids,
		map<string, map<double, vector<int>>> &code_ids, double th);
bool compactsc(std::vector<NSPgeometry::XYZ> sccrds, Subsite sub_new);
void rankformssites(std::vector<std::shared_ptr<Subsite>> &ssites,
		double thread, std::map<double, std::vector<int>> test_scores,
		std::vector<std::shared_ptr<Subsite>> test_ssites);
vector<string> maxscorecodes(std::map<std::string, std::map<double, std::vector<int>>> fcode_ids, string tot_code);
int unirandid(std::map<double, std::vector<int>> fci);
// return map<-bf_num, <tmpltmols>>
std::map<int, std::vector<OBMol>> readalltmplts(std::set<std::string> tmplt_names,
		OBConversion obconversion, OBMol targetmol, double perc);

// read tmplt_OBMol from name, and rank them by atm_cover_num. return
// <-atm_contact_num, vector<vector<id>>>
// because we can know atm_code and th here, we can directly know the minim_tmplt to be read.
vector<vector<int>> readandrank(std::vector<OBMol> &tmplt_ombs,
		std::vector<std::shared_ptr<SubstrAlignment>> &tmplt_ssas,
		std::vector<string> &tmplt_acode, std::set<std::string> tmplt_names,
		OBConversion ob, OBMol targetmol, double perc, int needs, double th);
// fetch the right ssa matches the code
std::vector<std::shared_ptr<SubstrAlignment>> ssaincode(string code, std::vector<std::shared_ptr<SubstrAlignment>> ssas);
// pick the top th local_ssites and push_back to ssites.
void pickssites(std::vector<TmpltSSAs> &ssas, TargetStruct &tgt, double weight,
		std::vector<std::shared_ptr<Subsite>> &ssites, int pick);
void topsimitmplts(std::vector<OBMol> &tmplt_mols, OBMol targetmol,
		std::vector<std::vector<std::string>> bfns, int needs, double perc);
}





#endif /* BUILDPOCKET_H_ */
