/*
 * linkerbuildercontrols.cpp
 *
 *  Created on: 2016年11月29日
 *      Author: hyliu
 */

#include <backbone/linkerbuildercontrols.h>
#include "dstl/randomengine.h"
#include "pdbstatistics/phipsidistr.h"

using namespace NSPproteinrep;
using namespace NSPdataio;

void LinkerBuilderControls::initcontrols(const std::string &filename) {
	ParameterSet & s=Controls::getparameterset();
	std::vector<std::string> intkeynames{"NCandidateLinkers","RandomSeed"};
	s.InitIntKeys(intkeynames);
	std::vector<std::string> stringkeynames{"SegmentsFile", "CoilPhiPsiFile","GlyPhiPsiFile","PreProPhiPsiFile",
	"ResultsFile","IdealGeometriesFile"};
	s.InitStringKeys(stringkeynames);

	std::vector<std::string> doublekeynames{};
	s.InitDoubleKeys(doublekeynames);

	std::vector<std::string> dveckeys{};
	s.InitDoubleVecKeys(dveckeys);

	std::vector<std::string> intveckeys{"LinkerLengthsAndCopies"};
	s.InitIntVecKeys(intveckeys);

	s.readparameters(filename);
	if(s.keydefined("RandomSeed")){
		int seed;
		s.getval("RandomSeed",&seed);
		NSPdstl::RandomEngine<>::getinstance().init((unsigned int)(seed));
	} else {
		NSPdstl::RandomEngine<>::getinstance().init();
	}
	LinkerBuilderControls & instance = LinkerBuilderControls::getinstance();
	int nc;
	s.getval("NCandidateLinkers", &nc);
	instance.ncandidates=nc;
	s.getval("SegmentsFile", &(instance.segmentsfilename));
	s.getval("ResultsFile", &(instance.outlinkerfilename));
	if(s.keydefined("IdealGeometriesFile")){
		s.getval("IdealGeometriesFile",&(instance.idealgeometriesfilename));
	} else {
		instance.idealgeometriesfilename="idealgeometries.dat";
	}
	if(s.keydefined("CoilPhiPsiFile")){
		std::string file;
		s.getval("CoilPhiPsiFile", &file);
		NSPpdbstatistics::PhiPsiDistr::coildistr(file);
	}
	if(s.keydefined("GlyPhiPsiFile")){
		std::string file;
		s.getval("GlyPhiPsiFile", &file);
		NSPpdbstatistics::PhiPsiDistr::glydistr(file);
	}
	if(s.keydefined("PreProPhiPsiFile")){
		std::string file;
		s.getval("PreProPhiPsiFile", &file);
		NSPpdbstatistics::PhiPsiDistr::glydistr(file);
	}
	std::vector<int> ivec;
	s.getval("LinkerLengthsAndCopies", &ivec);
	for(int i=0; i<ivec.size(); i+=4){
		std::vector<unsigned int> lc;
		for(int m=0; m<4; ++m){
			unsigned int t=ivec[i+m];
			lc.push_back(t);
		}
		instance.lengthsandcopies.push_back(lc);
	}
}



