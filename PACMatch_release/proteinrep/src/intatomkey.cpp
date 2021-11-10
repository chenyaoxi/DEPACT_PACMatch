/*
 * intaomkey.cpp
 *
 *  Created on: 2016年11月15日
 *      Author: hyliu
 */

#include "proteinrep/intatomkey.h"
#include <cassert>

const std::map<std::string,unsigned int> NSPproteinrep::atomcodes{
	{"N",0}, {"H",1}, {"CA",2},{"CB",3},{"CG",4},
	{"CG1",9},{"CG2",5},{"SG",6},{"OG",7},{"OG1",8},
	{"CD1",10},{"CD2",11},{"ND1",12},{"OD1",13},{"OD2",14},
	{"CE",15},{"CE1",16},{"CE2",17},{"OE1",18},{"OE2",19},{"NE1",20},{"NE2",21},
	{"CZ",22},{"CD",23}, {"CE3",24},{"CZ2",25},{"CZ3",26},{"CH2",27},{"NH1",28},
	{"NH2",29},{"OH",30},{"NE",31},{"NZ",32},{"ND2",33},{"SD",34},
	{"C",78},{"O",79}};

const std::set<std::string> NSPproteinrep::mainchainatoms{
	"N","H","CA","C","O"
};
const std::set<std::string> NSPproteinrep::hbdonoratoms{"N"};
const std::set<std::string> NSPproteinrep::hbacceptoratoms{"O"};

const std::map<std::string, unsigned int> NSPproteinrep::residuenamecodes{
	{"ANY",0},{"ALA",1},{"CYS",2},{"ASP",3},{"GLU",4},{"PHE",5},{"GLY",6},{"HIS",7},
	{"ILE",8},{"LYS",9},{"LEU",10},{"MET",11},{"ASN",12},{"PRO",13},{"GLN",14},{"ARG",15},
	{"SER",16},{"THR",17},{"VAL",18},{"TRP",19},{"TYR",20},{"CYSH",21},{"ASPH",22},{"GLUH",23},
	{"HISH",24},{"LYSH",25},{"ARGH",26},{"CYSD",27},{"ASPD",28},{"GLUD",29},{"HISD",30},
	{"HISE",31},{"LYSD",32},{"ARGD",33},{"UMC",34},{"ZMC",35}};

const std::set<std::string> NSPproteinrep::positiveresidues{"HIS","LYS","HISH","ARG","ARGH"};
const std::set<std::string> NSPproteinrep::negativeresidues{"GLU","ASP","CYSD","GLUD","ASPD"};
const std::set<std::string> NSPproteinrep::polarresidues{"ASN","GLN",
		"CYS","CYSH","HISD","HISE","SER","THR","TYR","GLUH","ASPH","LYSD","ARGD"};
const std::set<std::string> NSPproteinrep::aromaticresidues{"TYR","PHE","TRP"};

const std::map<char, std::string> NSPproteinrep::residuecharnames{
	{'X',"UDA"},{'A',"ALA"},{'C',"CYS"},{'D',"ASP"},{'E',"GLU"},{'F',"PHE"},{'G',"GLY"},{'H',"HIS"},
	{'I',"ILE"},{'K',"LYS"},{'L',"LEU"},{'M',"MET"},{'N',"ASN"},{'P',"PRO"},{'Q',"GLN"},{'R',"ARG"},
	{'S',"SER"},{'T',"THR"},{'V',"VAL"},{'W',"TRP"},{'Y',"TYR"},{'U',"UMC"},{'Z',"ZMC"}};
